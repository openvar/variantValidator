import pytest
from unittest.mock import MagicMock, patch

from VariantValidator.modules.vvDatabase import Database


@pytest.fixture
def db():
    db = Database.__new__(Database)
    db.pool = None
    return db


def test_in_entries_none(db):
    db.query_with_fetchone = MagicMock(
        return_value=["none", "No data"]
    )

    result = db.in_entries("NM_1", "transcript_info")

    assert result["none"] == "none"


def test_in_entries_error(db):
    db.query_with_fetchone = MagicMock(
        return_value=["error", "broken"]
    )

    result = db.in_entries("NM_1", "transcript_info")

    assert result["error"] == "error"


def test_in_entries_success(db):
    db.query_with_fetchone = MagicMock(
        return_value=[
            "NM_1",
            "desc",
            "var",
            "1",
            "GENE",
            "UTA",
            "today",
            "false",
        ]
    )

    result = db.in_entries("NM_1", "transcript_info")

    assert result["accession"] == "NM_1"
    assert result["hgnc_symbol"] == "GENE"


@patch("VariantValidator.modules.vvDatabase.time.sleep")
def test_data_add_retry(mock_sleep, db):

    db.update_transcript_info_record = MagicMock()

    db.in_entries = MagicMock(
        side_effect=[
            {"none": "none"},
            {"none": "none"},
            {"accession": "NM_1"},
        ]
    )

    result = db.data_add(
        "NM_1",
        MagicMock()
    )

    assert result == {"accession": "NM_1"}
    assert mock_sleep.call_count == 2


def test_update_refseqgene_loci_insert(db):
    db.get_refseq_data_by_refseq_id = MagicMock(
        return_value=["none"]
    )
    db.insert_refseq_gene_data = MagicMock()

    db.update_refseqgene_loci(
        ["NG_1", "x", "y"]
    )

    db.insert_refseq_gene_data.assert_called_once()


def test_update_refseqgene_loci_update(db):
    db.get_refseq_data_by_refseq_id = MagicMock(
        return_value=["exists"]
    )
    db.update_refseq_gene_data = MagicMock()

    db.update_refseqgene_loci(
        ["NG_1", "x", "y"]
    )

    db.update_refseq_gene_data.assert_called_once()


def test_update_lrg_rs_lookup(db):
    db.get_refseq_id_from_lrg_id = MagicMock(
        return_value="none"
    )
    db.insert_refseq_gene_id_from_lrg_id = MagicMock()

    db.update_lrg_rs_lookup(["LRG_1"])

    db.insert_refseq_gene_id_from_lrg_id.assert_called_once()


def test_update_lrgt_rst(db):
    db.get_refseq_transcript_id_from_lrg_transcript_id = MagicMock(
        return_value="none"
    )
    db.insert_lrg_transcript_data = MagicMock()

    db.update_lrgt_rst(["LRG_1t1"])

    db.insert_lrg_transcript_data.assert_called_once()


def test_update_lrg_p_rs_p_lookup(db):
    db.get_refseq_protein_id_from_lrg_protein_id = MagicMock(
        return_value="none"
    )
    db.insert_lrg_protein_data = MagicMock()

    db.update_lrg_p_rs_p_lookup(
        "LRG_1p1",
        "NP_1"
    )

    db.insert_lrg_protein_data.assert_called_once()


def test_ref_type_assign():
    db = Database.__new__(Database)
    db.pool = None

    assert db.ref_type_assign("NC_000001.11") == ":g."
    assert db.ref_type_assign("NG_000001.1") == ":g."
    assert db.ref_type_assign("NM_000001.1") == ":c."
    assert db.ref_type_assign("NR_000001.1") == ":n."
    assert db.ref_type_assign("NP_000001.1") == ":p."


def test_ref_type_assign_lrg_transcript_c():
    db = Database.__new__(Database)
    db.pool = None

    db.get_refseq_transcript_id_from_lrg_transcript_id = MagicMock(
        return_value="NM_000001.1"
    )

    assert db.ref_type_assign("LRG_1t1") == ":c."


def test_ref_type_assign_lrg_transcript_n():
    db = Database.__new__(Database)
    db.pool = None

    db.get_refseq_transcript_id_from_lrg_transcript_id = MagicMock(
        return_value="NR_000001.1"
    )

    assert db.ref_type_assign("LRG_1t1") == ":n."


def test_ref_type_assign_lrg_protein():
    db = Database.__new__(Database)
    db.pool = None

    assert db.ref_type_assign("LRG_1_p") == ":p."


def test_ref_type_assign_unknown():
    db = Database.__new__(Database)
    db.pool = None

    with pytest.raises(Exception):
        db.ref_type_assign("XYZ")

def test_update_lrg_rs_lookup_existing(db):
    db.get_refseq_id_from_lrg_id = MagicMock(
        return_value="NG_123"
    )
    db.insert_refseq_gene_id_from_lrg_id = MagicMock()

    db.update_lrg_rs_lookup(["LRG_1"])

    db.insert_refseq_gene_id_from_lrg_id.assert_not_called()


def test_update_lrgt_rst_existing(db):
    db.get_refseq_transcript_id_from_lrg_transcript_id = MagicMock(
        return_value="NM_123"
    )
    db.insert_lrg_transcript_data = MagicMock()

    db.update_lrgt_rst(["LRG_1t1"])

    db.insert_lrg_transcript_data.assert_not_called()


def test_update_lrg_p_rs_p_lookup_existing(db):
    db.get_refseq_protein_id_from_lrg_protein_id = MagicMock(
        return_value="NP_123"
    )
    db.insert_lrg_protein_data = MagicMock()

    db.update_lrg_p_rs_p_lookup(
        "LRG_1p1",
        "NP_1"
    )

    db.insert_lrg_protein_data.assert_not_called()


def test_ref_type_assign_nt():
    db = Database.__new__(Database)
    db.pool = None

    assert db.ref_type_assign("NT_123456.1") == ":g."


def test_ref_type_assign_nw():
    db = Database.__new__(Database)
    db.pool = None

    assert db.ref_type_assign("NW_123456.1") == ":g."


def test_ref_type_assign_lrg_genomic():
    db = Database.__new__(Database)
    db.pool = None

    assert db.ref_type_assign("LRG_1") == ":g."

@patch("VariantValidator.modules.vvDatabase.utils.hgnc_rest")
def test_update_gene_stable_identifiers_insert(mock_hgnc_rest, db):
    docs = {
        "hgnc_id": "HGNC:123",
        "entrez_id": "456",
        "ensembl_gene_id": "ENSG000001",
        "omim_id": ["123456"],
        "ucsc_id": "uc001abc",
        "vega_id": "OTTHUMG000001",
        "ccds_id": ["CCDS1"],
        "location": "1p36",
        "name": "Gene name",
        "symbol": "GENE1",
    }

    mock_hgnc_rest.return_value = {
        "record": {
            "response": {
                "numFound": 1,
                "docs": [docs],
            }
        }
    }

    db.get_stable_gene_id_from_hgnc_id = MagicMock(
        return_value=["none", "No data"]
    )
    db.insert_gene_stable_ids = MagicMock()
    db.update_gene_stable_ids = MagicMock()

    result = db.update_gene_stable_identifiers("GENE1")

    db.insert_gene_stable_ids.assert_called_once()
    db.update_gene_stable_ids.assert_not_called()

    assert result["map_loc"] == "1p36"
    assert result["gene_name"] == "Gene name"
    assert result["prev"] is None


@patch("VariantValidator.modules.vvDatabase.utils.hgnc_rest")
def test_update_gene_stable_identifiers_update(mock_hgnc_rest, db):
    docs = {
        "hgnc_id": "HGNC:123",
        "symbol": "GENE1",
    }

    mock_hgnc_rest.return_value = {
        "record": {
            "response": {
                "numFound": 1,
                "docs": [docs],
            }
        }
    }

    db.get_stable_gene_id_from_hgnc_id = MagicMock(
        return_value=["HGNC:123", "exists"]
    )
    db.insert_gene_stable_ids = MagicMock()
    db.update_gene_stable_ids = MagicMock()

    db.update_gene_stable_identifiers("GENE1")

    db.update_gene_stable_ids.assert_called_once()
    db.insert_gene_stable_ids.assert_not_called()


@patch("VariantValidator.modules.vvDatabase.utils.hgnc_rest")
def test_update_gene_stable_identifiers_previous_symbol(mock_hgnc_rest, db):
    mock_hgnc_rest.side_effect = [
        {
            "record": {
                "response": {
                    "numFound": 0
                }
            }
        },
        {
            "error": "false",
            "record": {
                "response": {
                    "numFound": 1,
                    "docs": [
                        {
                            "symbol": "NEWGENE"
                        }
                    ]
                }
            }
        },
        {
            "record": {
                "response": {
                    "numFound": 1,
                    "docs": [
                        {
                            "hgnc_id": "HGNC:999",
                            "symbol": "NEWGENE"
                        }
                    ]
                }
            }
        }
    ]

    db.get_stable_gene_id_from_hgnc_id = MagicMock(
        return_value=["none", "No data"]
    )
    db.insert_gene_stable_ids = MagicMock()

    db.update_gene_stable_identifiers("OLDGENE")

    db.insert_gene_stable_ids.assert_called_once()


@patch("VariantValidator.modules.vvDatabase.utils.hgnc_rest")
def test_update_gene_stable_identifiers_symbol_not_found(mock_hgnc_rest, db):
    mock_hgnc_rest.side_effect = [
        {
            "record": {
                "response": {
                    "numFound": 0
                }
            }
        },
        {
            "error": "false",
            "record": {
                "response": {
                    "numFound": 0
                }
            }
        },
    ]

    db.insert_gene_stable_ids = MagicMock()
    db.update_gene_stable_ids = MagicMock()

    assert db.update_gene_stable_identifiers("NOT_A_GENE") is None

def test_query_with_fetchone_no_row():
    db = Database.__new__(Database)

    conn = MagicMock()
    cursor = MagicMock()

    cursor.fetchone.return_value = None

    db.get_conn = MagicMock(return_value=conn)
    db.get_cursor = MagicMock(return_value=cursor)

    result = Database.query_with_fetchone.__wrapped__(db, "NM_000001.1")

    assert result == ["none", "No data"]

    cursor.close.assert_called_once()
    conn.close.assert_called_once()


@patch("VariantValidator.modules.vvDatabase.utils.hgnc_rest")
def test_update_gene_stable_identifiers_unassigned(mock_hgnc_rest, db):
    mock_hgnc_rest.assert_not_called()

    assert db.update_gene_stable_identifiers("unassigned") is None


def test_update_refseqgene_loci_passes_second_argument(db):
    db.get_refseq_data_by_refseq_id = MagicMock(return_value=["exists"])
    db.update_refseq_gene_data = MagicMock()

    data = ["NG_000001.1", "x", "GRCh38"]

    db.update_refseqgene_loci(data)

    db.get_refseq_data_by_refseq_id.assert_called_once_with(
        "NG_000001.1",
        "GRCh38",
    )


def test_ref_type_assign_lrg_transcript_non_nm():
    db = Database.__new__(Database)

    db.get_refseq_transcript_id_from_lrg_transcript_id = MagicMock(
        return_value="NR_999999.1"
    )

    assert db.ref_type_assign("LRG_1t1") == ":n."


@patch("VariantValidator.modules.vvDatabase.utils.hgnc_rest")
def test_update_gene_stable_identifiers_missing_optional_fields(
    mock_hgnc_rest,
    db,
):
    mock_hgnc_rest.return_value = {
        "record": {
            "response": {
                "numFound": 1,
                "docs": [
                    {
                        "hgnc_id": "HGNC:1",
                        "symbol": "GENE1",
                    }
                ],
            }
        }
    }

    db.get_stable_gene_id_from_hgnc_id = MagicMock(
        return_value=["none", "No data"]
    )

    db.insert_gene_stable_ids = MagicMock()

    result = db.update_gene_stable_identifiers("GENE1")

    assert result == {
        "map_loc": None,
        "gene_name": None,
        "prev": None,
    }

    db.insert_gene_stable_ids.assert_called_once()


@patch("VariantValidator.modules.vvDatabase.utils.hgnc_rest")
def test_update_gene_stable_identifiers_existing_record_without_optional_fields(
    mock_hgnc_rest,
    db,
):
    mock_hgnc_rest.return_value = {
        "record": {
            "response": {
                "numFound": 1,
                "docs": [
                    {
                        "hgnc_id": "HGNC:1",
                        "symbol": "GENE1",
                    }
                ],
            }
        }
    }

    db.get_stable_gene_id_from_hgnc_id = MagicMock(
        return_value=["HGNC:1", "exists"]
    )

    db.update_gene_stable_ids = MagicMock()

    db.update_gene_stable_identifiers("GENE1")

    db.update_gene_stable_ids.assert_called_once()


def test_data_add_returns_after_first_success(db):
    db.update_transcript_info_record = MagicMock()

    db.in_entries = MagicMock(
        return_value={"accession": "NM_000001.1"}
    )

    result = db.data_add(
        "NM_000001.1",
        MagicMock(),
    )

    assert result["accession"] == "NM_000001.1"

    assert db.in_entries.call_count == 1

@patch("VariantValidator.modules.vvDatabase.utils.hgnc_rest")
def test_update_gene_stable_identifiers_missing_hgnc_id(
    mock_hgnc_rest,
    db,
):
    mock_hgnc_rest.return_value = {
        "record": {
            "response": {
                "numFound": 1,
                "docs": [
                    {
                        "symbol": "GENE1",
                    }
                ]
            }
        }
    }

    db.get_stable_gene_id_from_hgnc_id = MagicMock(
        return_value=["none", "No data"]
    )
    db.insert_gene_stable_ids = MagicMock()

    db.update_gene_stable_identifiers("GENE1")

    inserted = db.insert_gene_stable_ids.call_args.args[0]

    assert inserted["hgnc_id"] == ""
    assert inserted["entrez_id"] == ""
    assert inserted["ensembl_gene_id"] == ""
    assert inserted["ucsc_id"] == ""
    assert inserted["vega_id"] == ""


@patch("VariantValidator.modules.vvDatabase.time.sleep")
def test_data_add_exhausts_retry_limit(
    mock_sleep,
    db,
):
    db.update_transcript_info_record = MagicMock()

    db.in_entries = MagicMock(
        return_value={"none": "none"}
    )

    result = db.data_add(
        "NM_1",
        MagicMock(),
    )

    assert result == {"none": "none"}

    assert db.in_entries.call_count == 10
    assert mock_sleep.call_count == 9


def test_in_entries_unknown_table(db):
    assert db.in_entries(
        "NM_1",
        "not_a_table",
    ) == {}


@patch("VariantValidator.modules.vvDatabase.utils.hgnc_rest")
def test_update_gene_stable_identifiers_previous_symbol_with_optional_fields(
    mock_hgnc_rest,
    db,
):
    mock_hgnc_rest.side_effect = [
        {
            "record": {
                "response": {
                    "numFound": 0
                }
            }
        },
        {
            "error": "false",
            "record": {
                "response": {
                    "numFound": 1,
                    "docs": [
                        {
                            "symbol": "GENE2"
                        }
                    ]
                }
            }
        },
        {
            "record": {
                "response": {
                    "numFound": 1,
                    "docs": [
                        {
                            "hgnc_id": "HGNC:2",
                            "symbol": "GENE2",
                            "prev_symbol": ["OLD1"],
                            "location": "2q37",
                            "name": "Gene Two",
                        }
                    ]
                }
            }
        },
    ]

    db.get_stable_gene_id_from_hgnc_id = MagicMock(
        return_value=["none", "No data"]
    )

    db.insert_gene_stable_ids = MagicMock()

    result = db.update_gene_stable_identifiers("OLD1")

    assert result["prev"] == ["OLD1"]
    assert result["map_loc"] == "2q37"
    assert result["gene_name"] == "Gene Two"

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