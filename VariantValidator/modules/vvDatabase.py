from . import utils
from .utils import handleCursor
from . import vvDBInsert
import re
import vvhgvs.exceptions
import logging
import json

logger = logging.getLogger(__name__)


class Database(vvDBInsert.Mixin):
    """
    This class contains and handles the mysql connections for the variant validator database.

    It now uses mixins, and the order of inheritance is
    vvDBInit.Mixin
       v
    vvDBGet.Mixin
       v
    vvDBInsert.Mixin
       v
    vvDatabase
    """
    # from dbquery
    @handleCursor
    def query_with_fetchone(self, entry):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        # Expiry set to 12 months because we will from 2021 be rolling out 3-monthly database dumps of the validator db
        query = "SELECT refSeqID, description, transcriptVariant, currentVersion, hgncSymbol, utaSymbol, updated, " \
                "IF(updated < NOW() - INTERVAL 12 MONTH , 'true', 'false') FROM transcript_info WHERE " \
                "refSeqID = '%s'" % entry
        cursor.execute(query)
        row = cursor.fetchone()
        if row is None:
            row = ['none', 'No data']
            logger.debug("No data returned from query " + str(query))

        # close conn
        cursor.close()
        conn.close()
        return row

    # From data
    def data_add(self, accession, validator):
        """
        # Add accurate transcript descriptions to the database
        :param accession:
        :param validator:
        :return:
        """
        self.update_transcript_info_record(accession, validator)
        entry = self.in_entries(accession, 'transcript_info')
        return entry

    def in_entries(self, entry, table):
        """
        Retrieve transcript information
        :param entry:
        :param table:
        :return:
        """
        data = {}
        if table == 'transcript_info':
            row = self.query_with_fetchone(entry)
            if row[0] == 'error':
                data['error'] = row[0]
                data['description'] = row[1]
            elif row[0] == 'none':
                data['none'] = row[0]
                data['description'] = row[1]
            else:
                data['accession'] = row[0]
                data['description'] = row[1]
                data['variant'] = row[2]
                data['version'] = row[3]
                data['hgnc_symbol'] = row[4]
                data['uta_symbol'] = row[5]
                data['updated'] = row[6]
                data['expiry'] = row[7]
        return data

    def update_gene_stable_identifiers(self, symbol):
        # First perform a search against the input gene symbol or the symbol inferred from UTA
        initial = utils.hgnc_rest(path="/fetch/symbol/" + symbol)

        # Check for a record
        if str(initial['record']['response']['numFound']) == '0':
            # Search hgnc rest to see if symbol is out of date
            rest_data = utils.hgnc_rest(path="/search/prev_symbol/" + symbol)
            # If the name is correct no record will be found
            if rest_data['error'] == 'false' and int(rest_data['record']['response']['numFound']) != 0:
                symbol = rest_data['record']['response']['docs'][0]['symbol']
                initial = utils.hgnc_rest(path="/fetch/symbol/" + symbol)

        if symbol != 'unassigned' and int(initial['record']['response']['numFound']) != 0:
            docs = initial['record']['response']['docs'][0]
            hgnc_id = ''
            entrez_id = ''
            ensembl_gene_id = ''
            omim_id = json.dumps([])
            ucsc_id = ''
            vega_id = ''
            ccds_id = json.dumps([])

            if 'hgnc_id' in docs:
                hgnc_id = docs['hgnc_id']
            if 'entrez_id' in docs:
                entrez_id = docs['entrez_id']
            if 'ensembl_gene_id' in docs:
                ensembl_gene_id = docs['ensembl_gene_id']
            if 'omim_id' in docs:
                omim_id = json.dumps(docs['omim_id'])
            if 'ucsc_id' in docs:
                ucsc_id = docs['ucsc_id']
            if 'vega_id' in docs:
                vega_id = docs['vega_id']
            if 'ccds_id' in docs:
                ccds_id = json.dumps(docs['ccds_id'])

            gene_stable_ids = {
                "hgnc_id": hgnc_id,
                "entrez_id": entrez_id,
                "ensembl_gene_id": ensembl_gene_id,
                "omim_id": omim_id,
                "ucsc_id": ucsc_id,
                "vega_id": vega_id,
                "ccds_id": ccds_id,
                "hgnc_symbol": symbol

            }
            gene_id_info = self.get_stable_gene_id_from_hgnc_id(gene_stable_ids["hgnc_id"])
            if gene_id_info[1] != 'No data':
                self.update_gene_stable_ids(gene_stable_ids)
            else:
                self.insert_gene_stable_ids(gene_stable_ids)

            # Compile additional data
            try:
                prev = docs["prev_symbol"]
            except KeyError:
                prev = None
            try:
                maploc = docs["location"]
            except KeyError:
                maploc = None
            try:
                name = docs["name"]
            except KeyError:
                name = None

            hgnc_data = {"map_loc": maploc,
                         "gene_name": name,
                         "prev": prev
            }

            return hgnc_data

    def update_transcript_info_record(self, accession, validator, bypass_with_symbol=False, **kwargs):

        """
        Search Ensembl APIs for transcript_info data
        """

        if 'ENST' in accession:

            """
            Ensembl APIs do not cross-reference GRCh38 and 38 so 37 queries are needed
            They also interchangeably decide whether or not they accept version information
            Therefore, assume they do not and check using Accession.Version Python split 
            """
            enst_accession, enst_version = accession.split('.')
            try:
                genome_build = kwargs['genome_build']
            except KeyError:
                genome_build = None

            # Make requests
            if genome_build is None:
                raise utils.DatabaseConnectionError("Connection to Ensembl database requires specification of "
                                                    "a genome build (GRCh37 or GRCh38)")
            elif genome_build is 'GRCh37' or genome_build is 'GRCh38':
                ens_record = utils.ensembl_rest(id=enst_accession, endpoint="/lookup/id/", genome=genome_build)
                ens_json = ens_record['record']

            # Check version
            try:
                if enst_version in str(ens_json['version']):
                    version = accession
                    description = str(ens_json['display_name'])
                    genbank_symbol = description.split('-')[0]
                    if ens_json['is_canonical'] == 1:
                        select_tx = 'Ensembl'
                    else:
                        select_tx = False
                    ensemblgene_id = ens_json['Parent']
                    mapped_chr = ens_json['seq_region_name']
                    map_position = "chr%s:%s:%s" % (mapped_chr, ens_json['start'], ens_json['end'])

                    # Get CCDS ID
                    ccds_record = utils.ensembl_rest(id=enst_accession,
                                                     endpoint="/xrefs/id/",
                                                     genome=genome_build,
                                                     options='external_db=CCDS')
                    if len(ccds_record) == 1:
                        ccds_id = ccds_record[0]['display_id']
                    else:
                        ccds_id = None

                    # Get gene db_xref
                    hgnc_id = None
                    gene_name = None
                    gene_xrefs = utils.ensembl_rest(id=ensemblgene_id,
                                                    endpoint="/xrefs/id/",
                                                    genome=genome_build)
                    gene_xrefs_json = gene_xrefs['record']
                    for xref in gene_xrefs_json:
                        if xref['dbname'] == 'HGNC':
                            hgnc_id = xref['primary_id']
                            gene_name = xref['description']

                    # Get MANE status and Ensembl canonical status
                    mane_select = False
                    ensembl_select = False
                    mane_plus_clinical = False
                    if select_tx == "Ensembl":
                        ensembl_select = True
                    tark_record = utils.ensembl_tark(id=enst_accession + '.' + enst_version,
                                                     endpoint="/api/transcript/stable_id_with_version/")
                    tark_json = tark_record['record']

                    if 'mane_transcript_type' in tark_json['results'][0].keys():
                        if tark_json['results'][0]['mane_transcript_type'] == "MANE SELECT":
                            mane_select = True
                            select_tx = 'MANE'
                        else:
                            mane_plus_clinical = True

                    # Compile metadata dictionary
                    variant = {"db_xref": {"ensemblgene": ensemblgene_id,
                                           "ncbigene": None,
                                           "HGNC": hgnc_id,
                                           "CCDS": ccds_id,
                                           "select": select_tx},
                               "chromosome": mapped_chr,
                               "map": map_position,
                               # "gene_synonym": synonyms,
                               "note": gene_name,
                               "variant": description.split('-')[1],
                               "mane_select": mane_select,
                               "mane_plus_clinical": mane_plus_clinical,
                               "ensembl_select": ensembl_select,
                               "refseq_select": False}

                else:
                    warning = "Ensembl transcript %s is not identified in the Ensembl APIs" % accession
                    raise utils.DatabaseConnectionError(warning)
            except TypeError:
                connection_error = "Cannot retrieve data from Ensembl REST for record %s" % accession
                if bypass_with_symbol is not False:
                    try:
                        self.update_gene_stable_identifiers(bypass_with_symbol)
                    except Exception as e:
                        logger.debug("Except pass, %s", e)
                        logger.info("Unable to connect to genenames.org with symbol %s", bypass_with_symbol)
                        connection_error = "Cannot connect to genenames.org with symbol %s", bypass_with_symbol
                raise utils.DatabaseConnectionError(connection_error)
            except Exception as e:
                warning = "Ensembl transcript %s is not identified in the Ensembl APIs" % accession
                warning = "Ensembl transcript %s is not identified in the Ensembl APIs" % accession
                raise utils.DatabaseConnectionError(warning)
        else:
            """
            Search Entrez for corresponding record for the RefSeq ID
            """

            try:
                record = validator.entrez_efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            except IOError:
                connection_error = "Cannot retrieve data from NCBI Entrez for record %s" % accession
                if bypass_with_symbol is not False:
                    try:
                        self.update_gene_stable_identifiers(bypass_with_symbol)
                    except Exception as e:
                        logger.debug("Except pass, %s", e)
                        logger.info("Unable to connect to genenames.org with symbol %s", bypass_with_symbol)
                        connection_error = "Cannot connect to genenames.org with symbol %s", bypass_with_symbol
                raise utils.DatabaseConnectionError(connection_error)

            version = record.id
            description = record.description
            try:
                genbank_symbol = str(record.features[1].qualifiers['gene'][0])
            except KeyError:
                raise utils.DatabaseConnectionError("Gene information is not available in the RefSeq record. Record "
                                                    "potentially deprecated")
            try:
                # Genbank can be out-of-date so check this is not a historic record
                # First perform a search against the input gene symbol or the symbol inferred from UTA
                initial = utils.hgnc_rest(path="/fetch/symbol/" + genbank_symbol)
                # Check for a record
                if str(initial['record']['response']['numFound']) != '0':
                    genbank_symbol = genbank_symbol
                # No record found, is it a previous symbol?
                else:
                    # Look up current name
                    current = utils.hgnc_rest(path="/search/prev_symbol/" + genbank_symbol)
                    # Look for historic names
                    # If historic names = 0
                    if str(current['record']['response']['numFound']) == '0':
                        genbank_symbol = genbank_symbol
                    else:
                        genbank_symbol = current['record']['response']['docs'][0]['symbol']
            except Exception:
                pass

            if 'transcript variant' in description:
                tv = re.search(r'transcript variant \w+', description)
                tv = str(tv.group(0))
                tv = tv.replace('transcript variant', '')
                variant = tv.strip()
                variant = variant.upper()  # Some tv descriptions are a or A
            else:
                variant = '0'

            # Add tags
            # Currently used format transcriptVariant|MANE
            my_quals = record.features[0].qualifiers
            my_quals = json.loads(json.dumps(my_quals))
            my_tags = record.features[1].qualifiers
            my_tags = json.loads(json.dumps(my_tags))
            my_tags['variant'] = variant

            for ft in record.features:
                if 'db_xref' in ft.qualifiers:
                    for ccds_id_is in ft.qualifiers['db_xref']:
                        if 'CCDS:' in ccds_id_is:
                            my_tags['db_xref'] = my_tags['db_xref'] + [ccds_id_is]

            for keywd in record.annotations['keywords']:
                if 'Select' in keywd:
                    my_tags['db_xref'] = my_tags['db_xref'] + ['select:' + keywd.split(' ')[0]]
                else:
                    my_tags['db_xref'] = my_tags['db_xref'] + ['select:' + 'False']

            # Dict the db_xref
            db_xrefs_dict = {}
            db_xrefs = my_tags['db_xref']
            for xref in db_xrefs:
                if not 'HGNC' in xref:
                    tag, label = xref.split(':')
                    db_xrefs_dict[tag] = label
                else:
                    tag_label = xref.split(':')
                    db_xrefs_dict[tag_label[1]] = tag_label[2]

            my_tags['db_xref'] = db_xrefs_dict

            # merge dicts
            all_tags = {**my_quals, **my_tags}

            # Format dict
            all_tags_formatted = {}
            for key, val in all_tags.items():
                if key == 'gene' or key == 'mol_type' or key == 'organism':
                    continue
                if len(val) == 1:
                    val = val[0]
                    if ';' in val:
                        val = val.replace(' ', '')
                        val = val.split(';')
                all_tags_formatted[key] = val

            # Compile metadata dictionary
            all_tags_formatted["db_xref"]["ncbigene"] = all_tags_formatted["db_xref"].pop("GeneID")
            all_tags_formatted["db_xref"]["ensemblgene"] = None
            all_tags_formatted["refseq_select"] = False
            all_tags_formatted["mane_select"] = False
            all_tags_formatted["ensembl_select"] = False
            all_tags_formatted["mane_plus_clinical"] = False  # Placeholder, not seen in RefSeq records yet
            if all_tags_formatted["db_xref"]["select"] == "MANE":
                all_tags_formatted["mane_select"] = True
                all_tags_formatted["refseq_select"] = True  # Assumes no conflict between MANE Select and RefSeq Select
            if all_tags_formatted["db_xref"]["select"] == "RefSeq":
                all_tags_formatted["refseq_select"] = True
            try:
                all_tags_formatted["db_xref"]["HGNC"] = "HGNC:" + all_tags_formatted["db_xref"]["HGNC"]
            except KeyError:
                all_tags_formatted["db_xref"]["HGNC"] = None
            if all_tags_formatted["db_xref"]["select"] == "False":
                all_tags_formatted["db_xref"]["select"] = False

            variant = all_tags_formatted

        """
        Get information from UTA
        """
        if "test" not in kwargs.keys() or kwargs["test"] is not True:
            try:
                uta_info = validator.hdp.get_tx_identity_info(version)
            except vvhgvs.exceptions.HGVSDataNotAvailableError:
                version_ac_ver = version.split('.')
                version = version_ac_ver[0] + '.' + str(int(version_ac_ver[1]) - 1)
                try:
                    uta_info = validator.hdp.get_tx_identity_info(version)
                except vvhgvs.exceptions.HGVSDataNotAvailableError:
                    raise utils.DatabaseConnectionError("The requested transcript was not found in the VVTA database")

            uta_symbol = str(uta_info[6])
            if uta_symbol == '':
                uta_symbol = 'unassigned'
        else:
            uta_symbol = 'unassigned'

        """
        Update Gene stable IDs from HGNC rest, genenames.org and metadata fields
        """
        # Delete "gene_synonym" key
        try:
            del variant["gene_synonym"]
        except KeyError:
            pass
        variant["previous_symbol"] = None
        try:
            hgnc_data = self.update_gene_stable_identifiers(genbank_symbol)
        except Exception as e:
            logger.debug("Except pass, %s", e)
            logger.info("Unable to connect to HGNC with symbol %s", genbank_symbol)
        if hgnc_data is not None:
            variant["map"] = hgnc_data["map_loc"]
            variant["note"] = hgnc_data["gene_name"]
            variant["previous_symbol"] = hgnc_data["prev"]
        else:
            pass

        # Fill in missing keys
        try:
            variant["db_xref"]["hgnc"] = variant["db_xref"].pop("HGNC")
        except KeyError:
            variant["db_xref"]["hgnc"] = None
        try:
            variant["db_xref"].pop("MIM")
        except KeyError:
            pass
        try:
            variant["db_xref"]["CCDS"]
        except KeyError:
            variant["db_xref"]["CCDS"] = None
        try:
            variant.pop("previous_symbol")
        except KeyError:
            pass

        # print(json.dumps(variant, sort_keys=False, indent=4, separators=(',', ': ')))

        # Make into a json for storage
        variant = json.dumps(variant)

        """
        Insert/Update the transcript information
        """
        # Query information
        query_info = [version, description, variant, version, genbank_symbol, uta_symbol]
        table = 'transcript_info'

        # Update the transcript_info table
        returned_data = self.in_entries(version, table)
        # If the entry is not in the database add it
        if 'none' in returned_data:
            self.insert(version, query_info, table)
        # If the data in the entry has changed, update it
        else:
            self.update(version, query_info)

        return

    def update_refseqgene_loci(self, rsg_data):
        # First query the database
        entry_exists = self.get_refseq_data_by_refseq_id(rsg_data[0], rsg_data[2])
        if entry_exists[0] == 'none':
            self.insert_refseq_gene_data(rsg_data)
        else:
            self.update_refseq_gene_data(rsg_data)

    def update_lrg_rs_lookup(self, lrg_rs_lookup):
        # First query the database
        rsg_id = self.get_refseq_id_from_lrg_id(lrg_rs_lookup[0])
        if rsg_id == 'none':
            self.insert_refseq_gene_id_from_lrg_id(lrg_rs_lookup)

    def update_lrgt_rst(self, lrgtx_to_rst_id):
        # First query the database
        rst_id = self.get_refseq_transcript_id_from_lrg_transcript_id(lrgtx_to_rst_id[0])
        if rst_id == 'none':
            self.insert_lrg_transcript_data(lrgtx_to_rst_id)

    def update_lrg_p_rs_p_lookup(self, lrg_p, rs_p):
        # First query the database
        rsp_id = self.get_refseq_protein_id_from_lrg_protein_id(lrg_p)
        if rsp_id == 'none':
            self.insert_lrg_protein_data(lrg_p, rs_p)

    def ref_type_assign(self, accession):
        if 'NC_' in accession or 'NG_' in accession or 'NT_' in accession or 'NW_' in accession:
            ref_type = ':g.'
        elif accession.startswith('NM_'):
            ref_type = ':c.'
        elif accession.startswith('NR_'):
            ref_type = ':n.'
        elif accession.startswith('NP_'):
            ref_type = ':p.'
        elif accession.startswith('LRG_'):
            if 't' in accession:
                refseqtranscript_reference = self.get_refseq_transcript_id_from_lrg_transcript_id(accession)
                if refseqtranscript_reference.startswith('NM_'):
                    ref_type = ':c.'
                else:
                    ref_type = ':n.'
            elif '_p' in accession:
                ref_type = ':p.'
            else:
                ref_type = ':g.'
        else:
            # shouldn't reach this point
            raise Exception('Unable to recognise accession')
        return ref_type

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
