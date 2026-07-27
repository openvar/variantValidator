from . import utils
from .utils import handleCursor
from . import vvDBInsert

import json
import logging
import time

import vvhgvs.exceptions


logger = logging.getLogger(__name__)


class Database(vvDBInsert.Mixin):
    """
    This class contains and handles the MySQL connections for the VariantValidator database.

    It now uses mixins, and the order of inheritance is
    vvDBInit.Mixin
       v
    vvDBGet.Mixin
       v
    vvDBInsert.Mixin
       v
    vvDatabase
    """

    # From dbquery
    @handleCursor
    def query_with_fetchone(self, entry):
        # Connect and create cursor
        conn = self.get_conn()
        cursor = self.get_cursor(conn)

        # Expiry set to 12 months because we will from 2021 be rolling out
        # 3-monthly database dumps of the validator db
        query = (
            "SELECT refSeqID, description, transcriptVariant, currentVersion, "
            "hgncSymbol, utaSymbol, updated, "
            "IF(updated < NOW() - INTERVAL 12 MONTH, 'true', 'false') "
            "FROM transcript_info WHERE refSeqID = %s"
        )
        cursor.execute(query, (entry,))
        row = cursor.fetchone()

        if row is None:
            row = ['none', 'No data']
            logger.debug("No data returned from query %s", query)

        cursor.close()
        conn.close()

        return row

    # From data
    def data_add(self, accession, validator, genome_build=None):
        """
        Add accurate transcript descriptions to the database.

        :param accession:
        :param validator:
        :param genome_build:
        :return:
        """
        self.update_transcript_info_record(
            accession,
            validator,
            genome_build=genome_build
        )

        entry = self.in_entries(accession, 'transcript_info')
        i = 1

        while i in range(10):
            if 'none' not in entry:
                break

            i += 1
            time.sleep(2)
            entry = self.in_entries(accession, 'transcript_info')

        return entry

    def in_entries(self, entry, table):
        """
        Retrieve transcript information.

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
        # First perform a search against the input gene symbol or the symbol
        # inferred from UTA
        initial = utils.hgnc_rest(path="/fetch/symbol/" + symbol)

        # Check whether the symbol is out of date
        if initial['record']['response']['numFound'] == 0:
            rest_data = utils.hgnc_rest(path="/search/prev_symbol/" + symbol)

            if (
                rest_data['error'] == 'false'
                and rest_data['record']['response']['numFound'] != 0
            ):
                symbol = rest_data['record']['response']['docs'][0]['symbol']
                initial = utils.hgnc_rest(path="/fetch/symbol/" + symbol)

        if (
            symbol != 'unassigned'
            and initial['record']['response']['numFound'] != 0
        ):
            docs = initial['record']['response']['docs'][0]

            hgnc_id = docs.get('hgnc_id', '')
            entrez_id = docs.get('entrez_id', '')
            ensembl_gene_id = docs.get('ensembl_gene_id', '')
            omim_id = json.dumps(docs.get('omim_id', []))
            ucsc_id = docs.get('ucsc_id', '')
            vega_id = docs.get('vega_id', '')
            ccds_id = json.dumps(docs.get('ccds_id', []))

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

            gene_id_info = self.get_stable_gene_id_from_hgnc_id(hgnc_id)

            if gene_id_info[1] != 'No data':
                self.update_gene_stable_ids(gene_stable_ids)
            else:
                self.insert_gene_stable_ids(gene_stable_ids)

            return {
                "map_loc": docs.get("location"),
                "gene_name": docs.get("name"),
                "prev": docs.get("prev_symbol"),
                "hgnc_id": hgnc_id
            }

    def update_transcript_info_record(
            self,
            accession,
            validator,
            bypass_with_symbol=False,
            **kwargs
    ):
        """
        Search Ensembl APIs or Entrez for transcript_info data.
        """

        if accession.startswith("ENST"):
            """
            Ensembl APIs do not cross-reference GRCh38 and 37, so GRCh37
            queries are needed. They also interchangeably decide whether or
            not they accept version information. Therefore, assume they do
            not and check using Accession.Version Python split.
            """
            enst_accession, enst_version = accession.strip().split(".")
            genome_build = kwargs.get("genome_build")

            if genome_build is None:
                raise utils.DatabaseConnectionError(
                    "Connection to Ensembl database requires specification "
                    "of a genome build (GRCh37 or GRCh38)"
                )

            if genome_build in ("GRCh37", "GRCh38"):
                ens_record = utils.ensembl_rest(
                    id=enst_accession,
                    endpoint="/lookup/id/",
                    genome=genome_build,
                )
                ens_json = ens_record["record"]

            try:
                if enst_version == str(ens_json["version"]):
                    version = accession.strip()
                    description = ens_json["display_name"]
                    genbank_symbol = "-".join(description.split("-")[:-1])

                    if ens_json["is_canonical"] == 1:
                        select_tx = "Ensembl"
                    else:
                        select_tx = False

                    ensemblgene_id = ens_json["Parent"]
                    mapped_chr = ens_json["seq_region_name"]
                    map_position = (
                        f"chr{mapped_chr}:"
                        f"{ens_json['start']}:"
                        f"{ens_json['end']}"
                    )

                    # Get CCDS ID.
                    ccds_record = utils.ensembl_rest(
                        id=enst_accession,
                        endpoint="/xrefs/id/",
                        genome=genome_build,
                        options="external_db=CCDS",
                    )

                    if len(ccds_record) == 1:
                        ccds_id = ccds_record[0]["display_id"]
                    else:
                        ccds_id = None

                    # Get gene db_xref.
                    hgnc_id = None
                    gene_name = None

                    gene_xrefs = utils.ensembl_rest(
                        id=ensemblgene_id,
                        endpoint="/xrefs/id/",
                        genome=genome_build,
                    )
                    gene_xrefs_json = gene_xrefs["record"]

                    for xref in gene_xrefs_json:
                        if xref["dbname"] == "HGNC":
                            hgnc_id = xref["primary_id"]
                            gene_name = xref["description"]

                    # Get MANE status and Ensembl canonical status.
                    mane_select = False
                    ensembl_select = select_tx == "Ensembl"
                    mane_plus_clinical = False

                    tark_record = utils.ensembl_tark(
                        id=f"{enst_accession}.{enst_version}",
                        endpoint="/api/transcript/stable_id_with_version/",
                    )
                    tark_json = tark_record["record"]

                    try:
                        mane_type = tark_json["results"][0].get(
                            "mane_transcript_type"
                        )
                    except IndexError:
                        mane_type = None

                    if mane_type == "MANE SELECT":
                        mane_select = True
                        select_tx = "MANE"
                    elif mane_type is not None:
                        mane_plus_clinical = True

                    description_parts = description.split("-")
                    if len(description_parts) == 3:
                        description = (
                            f"{description_parts[1]}-{description_parts[2]}"
                        )

                    # Compile metadata dictionary.
                    variant = {
                        "db_xref": {
                            "ensemblgene": ensemblgene_id,
                            "ncbigene": None,
                            "HGNC": hgnc_id,
                            "CCDS": ccds_id,
                            "select": select_tx,
                        },
                        "chromosome": mapped_chr,
                        "map": map_position,
                        "note": gene_name,
                        "variant": description.split("-")[1],
                        "mane_select": mane_select,
                        "mane_plus_clinical": mane_plus_clinical,
                        "ensembl_select": ensembl_select,
                        "refseq_select": False,
                    }

                else:
                    logger.warning(
                        "Version Mismatch for %s version=%s and record version=%s",
                        enst_accession,
                        enst_version,
                        ens_json["version"],
                    )
                    raise utils.DatabaseConnectionError(
                        f"Ensembl transcript {accession} is not identified in "
                        "the Ensembl APIs"
                    )

            except TypeError:
                connection_error = (
                    f"Cannot retrieve data from Ensembl REST for record "
                    f"{accession}"
                )

                if bypass_with_symbol is not False:
                    try:
                        self.update_gene_stable_identifiers(
                            bypass_with_symbol
                        )
                    except Exception as e:
                        logger.debug("Except pass, %s", e)
                        logger.info(
                            "Unable to connect to genenames.org with symbol %s",
                            bypass_with_symbol,
                        )
                        connection_error = (
                            f"Cannot connect to genenames.org with symbol "
                            f"{bypass_with_symbol}"
                        )

                raise utils.DatabaseConnectionError(connection_error)

            except Exception:
                raise utils.DatabaseConnectionError(
                    f"Ensembl transcript {accession} is not identified in "
                    "the Ensembl APIs"
                )

        else:
            """
            Search Entrez for corresponding record for the RefSeq ID.
            """
            try:
                record = validator.entrez_efetch(
                    db="nucleotide",
                    id=accession,
                    rettype="gb",
                    retmode="text",
                )
            except IOError:
                connection_error = (
                    f"Cannot currently retrieve data from NCBI Entrez for "
                    f"record {accession}"
                )

                if bypass_with_symbol is not False:
                    try:
                        self.update_gene_stable_identifiers(
                            bypass_with_symbol
                        )
                    except Exception as e:
                        logger.debug("Except pass, %s", e)
                        logger.info(
                            "Unable to connect to genenames.org with symbol %s",
                            bypass_with_symbol,
                        )
                        connection_error = (
                            f"Cannot connect to genenames.org with symbol "
                            f"{bypass_with_symbol}"
                        )

                raise utils.DatabaseConnectionError(connection_error)

            version = record.id
            description = record.description

            try:
                genbank_symbol = record.features[1].qualifiers["gene"][0]
            except KeyError:
                raise utils.DatabaseConnectionError(
                    "Gene information is not available in the RefSeq record. "
                    "Record potentially deprecated"
                )

            try:
                # GenBank can be out of date, so check whether this is a
                # historic gene symbol.
                initial = utils.hgnc_rest(
                    path=f"/fetch/symbol/{genbank_symbol}"
                )

                if initial["record"]["response"]["numFound"] == 0:
                    current = utils.hgnc_rest(
                        path=f"/search/prev_symbol/{genbank_symbol}"
                    )

                    if current["record"]["response"]["numFound"] != 0:
                        genbank_symbol = (
                            current["record"]["response"]["docs"][0]["symbol"]
                        )
            except Exception:
                pass

            if "transcript variant" in description:
                transcript_variant = description.split(
                    "transcript variant ",
                    1,
                )[1].split()[0]
                variant = transcript_variant.upper()
            else:
                variant = "0"

            # Add tags.
            # Currently used format transcriptVariant|MANE.
            my_quals = json.loads(
                json.dumps(record.features[0].qualifiers)
            )
            my_tags = json.loads(
                json.dumps(record.features[1].qualifiers)
            )
            my_tags["variant"] = variant

            for feature in record.features:
                if "db_xref" in feature.qualifiers:
                    for db_xref in feature.qualifiers["db_xref"]:
                        if db_xref.startswith("CCDS:"):
                            my_tags["db_xref"].append(db_xref)

            for keyword in record.annotations["keywords"]:
                if "Select" in keyword:
                    my_tags["db_xref"].append(
                        f"select:{keyword.split(' ')[0]}"
                    )
                elif "Plus Clinical" in keyword:
                    my_tags["db_xref"].append(
                        f"select:{keyword.replace(' ', '_')}"
                    )
                else:
                    my_tags["db_xref"].append("select:False")

            # Convert db_xref to a dictionary.
            db_xrefs_dict = {}

            for xref in my_tags["db_xref"]:
                if xref.startswith("HGNC:"):
                    tag_label = xref.split(":")
                    db_xrefs_dict[tag_label[1]] = tag_label[2]
                else:
                    tag, label = xref.split(":")
                    db_xrefs_dict[tag] = label

            my_tags["db_xref"] = db_xrefs_dict

            # Merge dictionaries.
            all_tags = {**my_quals, **my_tags}

            # Format dictionary.
            all_tags_formatted = {}

            for key, val in all_tags.items():
                if key in ("gene", "mol_type", "organism"):
                    continue

                if len(val) == 1:
                    val = val[0]

                    if ";" in val:
                        val = val.replace(" ", "").split(";")

                all_tags_formatted[key] = val

            # Compile metadata dictionary.
            db_xref = all_tags_formatted["db_xref"]

            db_xref["ncbigene"] = db_xref.pop("GeneID")
            db_xref["ensemblgene"] = None

            all_tags_formatted["refseq_select"] = False
            all_tags_formatted["mane_select"] = False
            all_tags_formatted["ensembl_select"] = False
            all_tags_formatted["mane_plus_clinical"] = False

            if db_xref["select"] == "MANE":
                all_tags_formatted["mane_select"] = True
                # Assumes no conflict between MANE Select and RefSeq Select.
                all_tags_formatted["refseq_select"] = True

            if db_xref["select"] == "RefSeq":
                all_tags_formatted["refseq_select"] = True

            if db_xref["select"] == "MANE_Plus_Clinical":
                all_tags_formatted["mane_plus_clinical"] = True

            try:
                db_xref["HGNC"] = f"HGNC:{db_xref['HGNC']}"
            except KeyError:
                db_xref["HGNC"] = None

            if db_xref["select"] == "False":
                db_xref["select"] = False

            variant = all_tags_formatted

        """
        Get information from UTA.
        """
        if kwargs.get("test") is not True:
            try:
                uta_info = validator.hdp.get_tx_identity_info(
                    version
                )
            except vvhgvs.exceptions.HGVSDataNotAvailableError:
                version_ac_ver = version.split(".")
                version = (
                    f"{version_ac_ver[0]}."
                    f"{int(version_ac_ver[1]) - 1}"
                )

                try:
                    uta_info = validator.hdp.get_tx_identity_info(
                        version
                    )
                except vvhgvs.exceptions.HGVSDataNotAvailableError:
                    raise utils.DatabaseConnectionError(
                        "The requested transcript was not found in the "
                        "VVTA database"
                    )

            uta_symbol = uta_info[6]

            if not uta_symbol:
                uta_symbol = "unassigned"
        else:
            uta_symbol = "unassigned"

        """
        Update Gene stable IDs from HGNC REST, genenames.org and metadata fields.
        """
        variant.pop("gene_synonym", None)
        variant["previous_symbol"] = None

        try:
            hgnc_data = self.update_gene_stable_identifiers(
                genbank_symbol
            )
        except Exception as e:
            logger.debug("Except pass, %s", e)
            logger.info(
                "Unable to connect to HGNC with symbol %s",
                genbank_symbol,
            )
            hgnc_data = None

        if hgnc_data is not None:
            variant["map"] = hgnc_data["map_loc"]
            variant["note"] = hgnc_data["gene_name"]
            variant["previous_symbol"] = hgnc_data["prev"]

        # Fill in missing keys.
        variant["db_xref"]["hgnc"] = (
            variant["db_xref"].pop("HGNC", None)
        )
        variant["db_xref"].pop("MIM", None)
        variant["db_xref"].setdefault("CCDS", None)
        variant.pop("previous_symbol", None)

        # Make into JSON for storage.
        variant = json.dumps(variant)

        """
        Insert/update the transcript information.
        """
        query_info = [
            version,
            description,
            variant,
            version,
            genbank_symbol,
            uta_symbol,
        ]
        table = "transcript_info"

        returned_data = self.in_entries(
            version,
            table,
        )

        if "none" in returned_data:
            self.insert(
                version,
                query_info,
                table,
            )
        else:
            self.update(
                version,
                query_info,
            )

    def update_refseqgene_loci(self, rsg_data):
        entry_exists = self.get_refseq_data_by_refseq_id(
            rsg_data[0],
            rsg_data[2]
        )

        if entry_exists[0] == 'none':
            self.insert_refseq_gene_data(rsg_data)
        else:
            self.update_refseq_gene_data(rsg_data)

    def update_lrg_rs_lookup(self, lrg_rs_lookup):
        rsg_id = self.get_refseq_id_from_lrg_id(lrg_rs_lookup[0])

        if rsg_id == 'none':
            self.insert_refseq_gene_id_from_lrg_id(lrg_rs_lookup)

    def update_lrgt_rst(self, lrgtx_to_rst_id):
        rst_id = self.get_refseq_transcript_id_from_lrg_transcript_id(
            lrgtx_to_rst_id[0]
        )

        if rst_id == 'none':
            self.insert_lrg_transcript_data(lrgtx_to_rst_id)

    def update_lrg_p_rs_p_lookup(self, lrg_p, rs_p):
        rsp_id = self.get_refseq_protein_id_from_lrg_protein_id(lrg_p)

        if rsp_id == 'none':
            self.insert_lrg_protein_data(lrg_p, rs_p)

    def ref_type_assign(self, accession):
        if accession.startswith(('NC_', 'NG_', 'NT_', 'NW_')):
            return ':g.'

        if accession.startswith('NM_'):
            return ':c.'

        if accession.startswith('NR_'):
            return ':n.'

        if accession.startswith('NP_'):
            return ':p.'

        if accession.startswith('LRG_'):
            if 't' in accession:
                refseq_transcript = (
                    self.get_refseq_transcript_id_from_lrg_transcript_id(
                        accession
                    )
                )

                if refseq_transcript.startswith('NM_'):
                    return ':c.'

                return ':n.'

            if '_p' in accession:
                return ':p.'

            return ':g.'

        # Shouldn't reach this point
        raise Exception('Unable to recognise accession')


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
