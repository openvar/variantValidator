import pytest
from unittest.mock import MagicMock

from VariantValidator.modules.transcript_map_data import TranscriptMapData


@pytest.fixture
def mock_hdp():
    hdp = MagicMock()

    hdp.get_tx_mapping_options.return_value = [
        ["NM_000001.1", "NC_000001.11", "splign", True, 1],
        ["NM_000001.1", "NC_000002.12", "blat", False, -1],
    ]

    hdp.get_tx_exons.return_value = [
        {"alt_strand": 1, "exon": 1},
        {"alt_strand": 1, "exon": 2},
    ]

    return hdp


class TestTranscriptMapData:

    def test_mapping_options_fetch(self, mock_hdp):
        tmd = TranscriptMapData(mock_hdp)

        result = tmd.mapping_options("NM_000001.1")

        assert len(result) == 2
        mock_hdp.get_tx_mapping_options.assert_called_once()

    def test_mapping_options_cached(self, mock_hdp):
        tmd = TranscriptMapData(mock_hdp)

        tmd.mapping_options("NM_000001.1")
        tmd.mapping_options("NM_000001.1")

        mock_hdp.get_tx_mapping_options.assert_called_once()

    def test_mapping_options_requires_hdp(self):
        tmd = TranscriptMapData()

        with pytest.raises(ValueError):
            tmd.mapping_options("NM_000001.1")

    def test_map_strand_positive(self, mock_hdp):
        tmd = TranscriptMapData(mock_hdp)

        assert tmd.map_strand(
            "NM_000001.1",
            "NC_000001.11"
        ) == 1

    def test_map_strand_missing(self, mock_hdp):
        tmd = TranscriptMapData(mock_hdp)

        assert tmd.map_strand(
            "NM_000001.1",
            "NC_DOES_NOT_EXIST"
        ) is False

    def test_is_gapped_map_true(self, mock_hdp):
        tmd = TranscriptMapData(mock_hdp)

        assert tmd.is_gapped_map(
            "NM_000001.1",
            "NC_000001.11"
        ) is True

    def test_is_gapped_map_false(self, mock_hdp):
        tmd = TranscriptMapData(mock_hdp)

        assert tmd.is_gapped_map(
            "NM_000001.1",
            "NC_000002.12"
        ) is False

    def test_map_type(self, mock_hdp):
        tmd = TranscriptMapData(mock_hdp)

        assert tmd.map_type(
            "NM_000001.1",
            "NC_000001.11"
        ) == "splign"

    def test_map_type_missing(self, mock_hdp):
        tmd = TranscriptMapData(mock_hdp)

        assert tmd.map_type(
            "NM_000001.1",
            "UNKNOWN"
        ) is False

    def test_mapped_exons(self, mock_hdp):
        tmd = TranscriptMapData(mock_hdp)

        result = tmd.mapped_exons(
            "NM_000001.1",
            "NC_000001.11"
        )

        assert len(result) == 2
        mock_hdp.get_tx_exons.assert_called_once()

    def test_tx_exons_forward_strand(self, mock_hdp):
        tmd = TranscriptMapData(mock_hdp)

        result = tmd.tx_exons(
            "NM_000001.1",
            "NC_000001.11",
            "splign"
        )

        assert result[0]["exon"] == 1

    def test_tx_exons_reverse_strand(self, mock_hdp):
        mock_hdp.get_tx_exons.return_value = [
            {"alt_strand": -1, "exon": 1},
            {"alt_strand": -1, "exon": 2},
        ]

        tmd = TranscriptMapData(mock_hdp)

        result = tmd.tx_exons(
            "NM_000001.1",
            "NC_000001.11",
            "splign"
        )

        assert result[0]["exon"] == 2

    def test_tx_exons_type_error(self, mock_hdp):
        mock_hdp.get_tx_exons.return_value = "bad_data"

        tmd = TranscriptMapData(mock_hdp)

        result = tmd.tx_exons(
            "NM_000001.1",
            "NC_000001.11",
            "splign"
        )

        assert result == "error"

def test_mapping_options_cache_hit(mock_hdp):
    tmd = TranscriptMapData(mock_hdp)

    expected = [["cached"]]
    tmd.mapping_opts["NM_000001.1"] = expected

    assert tmd.mapping_options("NM_000001.1") == expected
    mock_hdp.get_tx_mapping_options.assert_not_called()


def test_map_strand_cache_hit(mock_hdp):
    tmd = TranscriptMapData()

    tmd.mapped_strands = {
        "NM_000001.1": {
            "NC_000001.11": 1
        }
    }

    assert tmd.map_strand(
        "NM_000001.1",
        "NC_000001.11"
    ) == 1


def test_map_strand_cache_miss(mock_hdp):
    tmd = TranscriptMapData()

    tmd.mapped_strands = {
        "NM_000001.1": {}
    }

    assert tmd.map_strand(
        "NM_000001.1",
        "NC_000001.11"
    ) is False


def test_is_gapped_map_cache_true():
    tmd = TranscriptMapData()

    tmd.gap_status = {
        "NM_000001.1": {
            "NC_000001.11": True
        }
    }

    assert tmd.is_gapped_map(
        "NM_000001.1",
        "NC_000001.11"
    ) is True


def test_is_gapped_map_cache_false():
    tmd = TranscriptMapData()

    tmd.gap_status = {
        "NM_000001.1": {}
    }

    assert tmd.is_gapped_map(
        "NM_000001.1",
        "NC_000001.11"
    ) is False


def test_map_type_cache_hit():
    tmd = TranscriptMapData()

    tmd.mapping_types = {
        "NM_000001.1": {
            "NC_000001.11": "splign"
        }
    }

    assert tmd.map_type(
        "NM_000001.1",
        "NC_000001.11"
    ) == "splign"


def test_map_type_cache_missing():
    tmd = TranscriptMapData()

    tmd.mapping_types = {
        "NM_000001.1": {}
    }

    assert tmd.map_type(
        "NM_000001.1",
        "NC_000001.11"
    ) == []


def test_mapped_exons_requires_hdp():
    tmd = TranscriptMapData()

    with pytest.raises(ValueError):
        tmd.mapped_exons(
            "NM_000001.1",
            "NC_000001.11"
        )
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
