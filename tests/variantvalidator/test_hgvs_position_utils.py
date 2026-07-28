import pytest
import vvhgvs.parser
from vvhgvs.enums import Datum

from VariantValidator.modules import hgvs_position_utils


class TestHgvsPositionUtils:

    @classmethod
    def setup_class(cls):
        cls.hp = vvhgvs.parser.Parser()

    # ------------------------------------------------------------------
    # Object validation
    # ------------------------------------------------------------------

    def test_is_object(self):
        variant = self.hp.parse_hgvs_variant(
            'NM_000088.4:c.100A>G'
        )
        assert hgvs_position_utils.is_object(variant) is None

    def test_is_object_rejects_string(self):
        with pytest.raises(
                hgvs_position_utils.HgvsPositionException,
                match='is not a parsed hgvs object'
        ):
            hgvs_position_utils.is_object(
                'NM_000088.4:c.100A>G'
            )

    # ------------------------------------------------------------------
    # Intronic positions
    # ------------------------------------------------------------------

    @pytest.mark.parametrize(
        'description,start_intronic,end_intronic',
        [
            ('NM_000088.4:c.100A>G', False, False),
            ('NM_000088.4:c.100+1A>G', True, True),
            ('NM_000088.4:c.100+1_101del', True, False),
            ('NM_000088.4:c.100_101-1del', False, True),
            ('NM_000088.4:c.100+1_101-1del', True, True),
        ]
    )
    def test_intronic_positions(
            self,
            description,
            start_intronic,
            end_intronic
    ):
        variant = self.hp.parse_hgvs_variant(description)

        assert (
            hgvs_position_utils.start_position_is_intronic(variant)
            is start_intronic
        )
        assert (
            hgvs_position_utils.end_position_is_intronic(variant)
            is end_intronic
        )
        assert (
            hgvs_position_utils.both_positions_are_intronic(variant)
            is (start_intronic and end_intronic)
        )
        assert (
            hgvs_position_utils.either_position_is_intronic(variant)
            is (start_intronic or end_intronic)
        )
        assert (
            hgvs_position_utils.only_one_position_is_intronic(variant)
            is (start_intronic != end_intronic)
        )

    # ------------------------------------------------------------------
    # Positive offsets
    # ------------------------------------------------------------------

    @pytest.mark.parametrize(
        'description,start_positive,end_positive',
        [
            ('NM_000088.4:c.100A>G', False, False),
            ('NM_000088.4:c.100+1A>G', True, True),
            ('NM_000088.4:c.100+1_101del', True, False),
            ('NM_000088.4:c.100_101+1del', False, True),
            ('NM_000088.4:c.100+1_101+2del', True, True),
        ]
    )
    def test_positive_offsets(
            self,
            description,
            start_positive,
            end_positive
    ):
        variant = self.hp.parse_hgvs_variant(description)

        assert (
            hgvs_position_utils.start_offset_is_positive(variant)
            is start_positive
        )
        assert (
            hgvs_position_utils.end_offset_is_positive(variant)
            is end_positive
        )
        assert (
            hgvs_position_utils.both_offsets_are_positive(variant)
            is (start_positive and end_positive)
        )
        assert (
            hgvs_position_utils.either_offset_is_positive(variant)
            is (start_positive or end_positive)
        )

    # ------------------------------------------------------------------
    # Negative offsets
    # ------------------------------------------------------------------

    @pytest.mark.parametrize(
        'description,start_negative,end_negative',
        [
            ('NM_000088.4:c.100A>G', False, False),
            ('NM_000088.4:c.100-1A>G', True, True),
            ('NM_000088.4:c.100-1_101del', True, False),
            ('NM_000088.4:c.100_101-1del', False, True),
            ('NM_000088.4:c.100-1_101-2del', True, True),
        ]
    )
    def test_negative_offsets(
            self,
            description,
            start_negative,
            end_negative
    ):
        variant = self.hp.parse_hgvs_variant(description)

        assert (
            hgvs_position_utils.start_offset_is_negative(variant)
            is start_negative
        )
        assert (
            hgvs_position_utils.end_offset_is_negative(variant)
            is end_negative
        )
        assert (
            hgvs_position_utils.both_offsets_are_negative(variant)
            is (start_negative and end_negative)
        )
        assert (
            hgvs_position_utils.either_offset_is_negative(variant)
            is (start_negative or end_negative)
        )

    # ------------------------------------------------------------------
    # 5-prime UTR
    # ------------------------------------------------------------------

    @pytest.mark.parametrize(
        'description,start_utr,end_utr',
        [
            ('NM_000088.4:c.100A>G', False, False),
            ('NM_000088.4:c.-10A>G', True, True),
            ('NM_000088.4:c.-10_10del', True, False),
            ('NM_000088.4:c.-10_-1del', True, True),
        ]
    )
    def test_5_prime_utr(
            self,
            description,
            start_utr,
            end_utr
    ):
        variant = self.hp.parse_hgvs_variant(description)

        assert (
            hgvs_position_utils.start_is_5_prime_utr(variant)
            is start_utr
        )
        assert (
            hgvs_position_utils.end_is_5_prime_utr(variant)
            is end_utr
        )

    # ------------------------------------------------------------------
    # 3-prime UTR inspection
    # ------------------------------------------------------------------

    @pytest.mark.parametrize(
        'description,start_utr,end_utr',
        [
            ('NM_000088.4:c.100A>G', False, False),
            ('NM_000088.4:c.*10A>G', True, True),
            ('NM_000088.4:c.100_*10del', False, True),
            ('NM_000088.4:c.*10_*20del', True, True),
        ]
    )
    def test_3_prime_utr(
            self,
            description,
            start_utr,
            end_utr
    ):
        variant = self.hp.parse_hgvs_variant(description)

        assert (
            hgvs_position_utils.start_is_3_prime_utr(variant)
            is start_utr
        )
        assert (
            hgvs_position_utils.end_is_3_prime_utr(variant)
            is end_utr
        )

    # ------------------------------------------------------------------
    # 3-prime UTR modification
    # ------------------------------------------------------------------

    def test_set_start_as_3_prime_utr(self):
        variant = self.hp.parse_hgvs_variant(
            'NM_000088.4:c.100_200del'
        )

        hgvs_position_utils.set_start_as_3_prime_utr(variant)

        assert variant.posedit.pos.start.datum == Datum.CDS_END
        assert variant.posedit.pos.end.datum != Datum.CDS_END

    def test_set_end_as_3_prime_utr(self):
        variant = self.hp.parse_hgvs_variant(
            'NM_000088.4:c.100_200del'
        )

        hgvs_position_utils.set_end_as_3_prime_utr(variant)

        assert variant.posedit.pos.start.datum != Datum.CDS_END
        assert variant.posedit.pos.end.datum == Datum.CDS_END

    def test_set_both_positions_as_3_prime_utr(self):
        variant = self.hp.parse_hgvs_variant(
            'NM_000088.4:c.100_200del'
        )

        hgvs_position_utils.set_both_positions_as_3_prime_utr(
            variant
        )

        assert variant.posedit.pos.start.datum == Datum.CDS_END
        assert variant.posedit.pos.end.datum == Datum.CDS_END

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
