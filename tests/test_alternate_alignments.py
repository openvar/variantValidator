from VariantValidator import Validator
from unittest import TestCase
import vvhgvs
from unittest.mock import MagicMock, patch
import pytest

from VariantValidator.modules.format_converters import (
    map_alt_intron_to_primary,
    AltPrimaryIntronError,
    AltPrimaryMappingError,
)


class TestAltAlignments(TestCase):

    @classmethod
    def setup_class(cls):
        cls.vv = Validator()
        cls.vv.testing = True

    def test_sub(self):
        variant = 'NW_011332691.1(NM_012234.6):c.335+1G>C'
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_012234.6:c.335+1G>C" in list(results.keys())
        assert results["NM_012234.6:c.335+1G>C"][
                   "genome_context_intronic_sequence"] == "NW_011332691.1(NM_012234.6):c.335+1G>C"
        assert results["NM_012234.6:c.335+1G>C"]["alt_genomic_loci"][0]["grch38"][
            "hgvs_genomic_description"] == "NW_011332691.1:g.53828C>G"

    def test_del(self):
        variant = 'NW_011332691.1(NM_012234.6):c.335+1del'
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_012234.6:c.335+1del" in list(results.keys())
        assert "NW_011332691.1(NM_012234.6):c.335+1del" in results[
            "NM_012234.6:c.335+1del"]["genome_context_intronic_sequence"]
        assert "NW_011332691.1:g.53828del" in results["NM_012234.6:c.335+1del"]["alt_genomic_loci"][0]["grch38"][
            "hgvs_genomic_description"]

    def test_delins(self):
        variant = "NW_011332691.1(NM_012234.6):c.335+1delinsAT"
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_012234.6:c.335+1delinsAT" in list(results.keys())
        assert "NW_011332691.1(NM_012234.6):c.335+1delinsAT" in results[
            "NM_012234.6:c.335+1delinsAT"]["genome_context_intronic_sequence"]
        assert "NW_011332691.1:g.53828delinsAT" in results["NM_012234.6:c.335+1delinsAT"]["alt_genomic_loci"][
            0]["grch38"]["hgvs_genomic_description"]

    def test_inv(self):
        variant = "NW_011332691.1(NM_012234.6):c.335+1_335+6inv"
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_012234.6:c.335+1_335+6inv" in list(results.keys())
        assert "NW_011332691.1(NM_012234.6):c.335+1_335+6inv" in results[
            "NM_012234.6:c.335+1_335+6inv"]["genome_context_intronic_sequence"]
        assert "NW_011332691.1:g.53823_53828inv" in results["NM_012234.6:c.335+1_335+6inv"]["alt_genomic_loci"][0][
            "grch38"]["hgvs_genomic_description"]

    def test_dup(self):
        variant = "NW_011332691.1(NM_012234.6):c.335+1dup"
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_012234.6:c.335+1dup" in list(results.keys())
        assert "NW_011332691.1(NM_012234.6):c.335+1dup" in results[
            "NM_012234.6:c.335+1dup"]["genome_context_intronic_sequence"]
        assert "NW_011332691.1:g.53828dup" in results["NM_012234.6:c.335+1dup"]["alt_genomic_loci"][0]["grch38"][
            "hgvs_genomic_description"]

    def test_identity(self):
        variant = "NW_011332691.1(NM_012234.6):c.335+1_335+6="
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert "NM_012234.6:c.335+1_335+6=" in list(results.keys())
        assert "NW_011332691.1(NM_012234.6):c.335+1_335+6=" in results[
            "NM_012234.6:c.335+1_335+6="]["genome_context_intronic_sequence"]
        assert "NW_011332691.1:g.53823_53828=" in results["NM_012234.6:c.335+1_335+6="]["alt_genomic_loci"][0]["grch38"][
            "hgvs_genomic_description"]

    def test_bad_positive_intron(self):
        variant = "NW_011332691.1(NM_012234.6):c.330+1_330+6="
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert results["validation_warning_1"]["validation_warnings"] == [
            'ExonBoundaryError: Position c.330+1 does not correspond with an exon boundary for transcript reference '
            'NM_012234.6 and genomic reference NW_011332691.1']

    def test_bad_negative_intron(self):
        variant = "NW_011332691.1(NM_012234.6):c.338-6_338-2="
        results = self.vv.validate(variant, 'GRCh38', 'all').format_as_dict(test=True)
        print(results)
        assert results["validation_warning_1"]["validation_warnings"] ==  [
            'ExonBoundaryError: Position c.338-6 does not correspond with an exon boundary for transcript reference '
            'NM_012234.6 and genomic reference NW_011332691.1']

    def test_nt_grch37_negative_intron(self):
        variant = "NT_167248.1(NM_032470.3):c.215-2_215-1del"
        results = self.vv.validate(
            variant, 'GRCh37', 'all'
        ).format_as_dict(test=True)
        print(results)

        assert "NM_032470.3:c.215-2_215-1del" in list(results.keys())
        assert results["NM_032470.3:c.215-2_215-1del"][
            "genome_context_intronic_sequence"
        ] == "NC_000006.11(NM_032470.3):c.215-2_215-1del"
        assert results["NM_032470.3:c.215-2_215-1del"][
            "hgvs_refseqgene_variant"
        ] == "NG_008337.2:g.69657_69658del"

    def test_nt_grch38_negative_intron(self):
        variant = "NT_167248.2(NM_032470.3):c.215-2_215-1dup"
        results = self.vv.validate(
            variant, 'GRCh38', 'all'
        ).format_as_dict(test=True)
        print(results)

        assert "NM_032470.3:c.215-2_215-1dup" in list(results.keys())
        assert results["NM_032470.3:c.215-2_215-1dup"][
            "genome_context_intronic_sequence"
        ] == "NC_000006.12(NM_032470.3):c.215-2_215-1dup"

    def test_nt_grch37_negative_intron_on_38(self):
        variant = "NT_167248.1(NM_032470.3):c.215-2_215-1del"
        results = self.vv.validate(
            variant, 'GRCh38', 'all'
        ).format_as_dict(test=True)
        print(results)

        assert "NM_032470.3:c.215-2_215-1del" in list(results.keys())
        assert results["NM_032470.3:c.215-2_215-1del"][
            "genome_context_intronic_sequence"
        ] == "NC_000006.12(NM_032470.3):c.215-2_215-1del"
        assert results["NM_032470.3:c.215-2_215-1del"][
            "hgvs_refseqgene_variant"
        ] == "NG_008337.2:g.69657_69658del"

    def test_nt_grch38_negative_intron_on_37(self):
        variant = "NT_167248.2(NM_032470.3):c.215-2_215-1dup"
        results = self.vv.validate(
            variant, 'GRCh37', 'all'
        ).format_as_dict(test=True)
        print(results)

        assert "NM_032470.3:c.215-2_215-1dup" in list(results.keys())
        assert results["NM_032470.3:c.215-2_215-1dup"][
            "genome_context_intronic_sequence"
        ] == "NC_000006.11(NM_032470.3):c.215-2_215-1dup"


@pytest.fixture
def alt_intron_fixture():
    variant = MagicMock()
    validator = MagicMock()

    hgvs_transcript = MagicMock()
    hgvs_transcript.ac = "NM_000001.1"
    hgvs_transcript.type = "n"

    hgvs_transcript.posedit.pos.start.base = 100
    hgvs_transcript.posedit.pos.start.offset = 1
    hgvs_transcript.posedit.pos.end.base = 100
    hgvs_transcript.posedit.pos.end.offset = 1

    hgvs_transcript.posedit.edit.type = "sub"
    hgvs_transcript.posedit.edit.ref = "A"
    hgvs_transcript.posedit.edit.alt = "G"

    variant.quibble = hgvs_transcript
    variant.genomic_context_ac = "NW_000001.1"
    variant.primary_assembly = "GRCh37"

    # Prevent construction of real AssemblyMapper instances.
    variant.no_norm_evm = MagicMock()
    variant.evm = MagicMock()

    variant.hn = MagicMock()
    variant.reverse_normalizer = MagicMock()

    variant.warnings = []
    variant.coding = None
    variant.protein = ""
    variant.gene_symbol = ""
    variant.output_type_flag = None

    variant.genome_context_intronic_sequence = None
    variant.refseqgene_context_intronic_sequence = None

    validator.alt_aln_method = "splign"
    validator.vm = MagicMock()
    validator.hdp = MagicMock()

    return variant, validator


def make_exons(
        alt_ac,
        *,
        intron=True,
        strand=1,
        boundary=100,
):
    """
    Return minimal UTA-like exon records for an n.100+1 boundary.

    intron=True:
        genomic gap exists between the two exons.

    intron=False:
        genomic exons are directly adjacent.
    """
    second_start = 1200 if intron else 1100

    return [
        {
            "tx_ac": "NM_000001.1",
            "alt_ac": alt_ac,
            "alt_aln_method": "splign",
            "alt_strand": strand,
            "ord": 0,
            "tx_start_i": 0,
            "tx_end_i": boundary,
            "alt_start_i": 1000,
            "alt_end_i": 1100,
            "cigar": "100=",
        },
        {
            "tx_ac": "NM_000001.1",
            "alt_ac": alt_ac,
            "alt_aln_method": "splign",
            "alt_strand": strand,
            "ord": 1,
            "tx_start_i": boundary,
            "tx_end_i": 200,
            "alt_start_i": second_start,
            "alt_end_i": second_start + 100,
            "cigar": "100=",
        },
    ]


def make_hgvs_genomic(
        ac="NW_000001.1",
        ref="A",
        alt="G",
):
    hgvs_genomic = MagicMock()
    hgvs_genomic.ac = ac
    hgvs_genomic.posedit.edit.ref = ref
    hgvs_genomic.posedit.edit.alt = alt
    return hgvs_genomic


def make_remapped_transcript():
    transcript = MagicMock()
    transcript.ac = "NM_000001.1"
    transcript.type = "n"
    transcript.posedit.edit.type = "sub"
    transcript.posedit.edit.ref = "A"
    transcript.posedit.edit.alt = "G"
    return transcript


class TestAltIntronMocked:

    def test_non_alt_genomic_context_returns_false(
            self, alt_intron_fixture
    ):
        variant, validator = alt_intron_fixture
        variant.genomic_context_ac = "NC_000001.11"

        result = map_alt_intron_to_primary(
            variant, validator
        )

        assert result is False
        validator.hdp.get_tx_exons.assert_not_called()


    def test_none_genomic_context_returns_false(
            self, alt_intron_fixture
    ):
        variant, validator = alt_intron_fixture
        variant.genomic_context_ac = None

        result = map_alt_intron_to_primary(
            variant, validator
        )

        assert result is False
        validator.hdp.get_tx_exons.assert_not_called()


    def test_non_intronic_returns_false(
            self, alt_intron_fixture
    ):
        variant, validator = alt_intron_fixture

        with patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.either_position_is_intronic",
            return_value=False,
        ):
            result = map_alt_intron_to_primary(
                variant, validator
            )

        assert result is False
        validator.hdp.get_tx_exons.assert_not_called()


    def test_submitted_alignment_data_unavailable(
            self, alt_intron_fixture
    ):
        variant, validator = alt_intron_fixture

        validator.hdp.get_tx_exons.return_value = []

        with patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.either_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.start_position_is_intronic",
            return_value=True,
        ):
            with pytest.raises(
                AltPrimaryIntronError,
                match="AlignmentDataError",
            ):
                map_alt_intron_to_primary(
                    variant, validator
                )

        assert variant.coding is variant.quibble


    def test_submitted_get_tx_exons_error(
            self, alt_intron_fixture
    ):
        variant, validator = alt_intron_fixture

        validator.hdp.get_tx_exons.side_effect = (
            vvhgvs.exceptions.HGVSError("test")
        )

        with patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.either_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.start_position_is_intronic",
            return_value=True,
        ):
            with pytest.raises(
                AltPrimaryIntronError,
                match="AlignmentDataError",
            ):
                map_alt_intron_to_primary(
                    variant, validator
                )


    def test_submitted_bad_boundary_raises_exon_boundary_error(
            self, alt_intron_fixture
    ):
        variant, validator = alt_intron_fixture

        validator.hdp.get_tx_exons.return_value = make_exons(
            "NW_000001.1",
            boundary=90,
        )

        with patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.either_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.start_position_is_intronic",
            return_value=True,
        ):
            with pytest.raises(
                AltPrimaryIntronError,
                match="ExonBoundaryError",
            ):
                map_alt_intron_to_primary(
                    variant, validator
                )

        assert len(variant.warnings) == 1
        assert "ExonBoundaryError" in variant.warnings[0]


    def test_current_primary_contains_intron(
            self, alt_intron_fixture
    ):
        variant, validator = alt_intron_fixture

        validator.hdp.get_tx_mapping_options.return_value = [
            [
                "NM_000001.1",
                "NC_000001.10",
                "splign",
            ],
        ]

        validator.hdp.get_tx_exons.side_effect = [
            make_exons("NW_000001.1"),
            make_exons("NC_000001.10"),
        ]

        expected_genomic = make_hgvs_genomic(
            "NC_000001.10"
        )

        validator.vm.t_to_g.return_value = expected_genomic

        with patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.either_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.start_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "seq_data.supported_for_mapping",
            side_effect=lambda ac, assembly:
                ac == "NC_000001.10"
                and assembly == "GRCh37",
        ):
            result = map_alt_intron_to_primary(
                variant, validator
            )

        assert result is False
        assert variant.hgvs_genomic is expected_genomic

        validator.vm.t_to_g.assert_called_once_with(
            variant.quibble,
            "NC_000001.10",
            alt_aln_method="splign",
        )


    def test_switches_to_other_primary_assembly(
            self, alt_intron_fixture
    ):
        variant, validator = alt_intron_fixture

        validator.hdp.get_tx_mapping_options.return_value = [
            [
                "NM_000001.1",
                "NC_000001.10",
                "splign",
            ],
            [
                "NM_000001.1",
                "NC_000001.11",
                "splign",
            ],
        ]

        validator.hdp.get_tx_exons.side_effect = [
            make_exons("NW_000001.1"),
            make_exons(
                "NC_000001.10",
                intron=False,
            ),
            make_exons("NC_000001.11"),
        ]

        expected_genomic = make_hgvs_genomic(
            "NC_000001.11"
        )

        validator.vm.t_to_g.return_value = expected_genomic

        with patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.either_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.start_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "seq_data.supported_for_mapping",
            side_effect=lambda ac, assembly: (
                ac == "NC_000001.10"
                and assembly == "GRCh37"
            ) or (
                ac == "NC_000001.11"
                and assembly == "GRCh38"
            ),
        ):
            result = map_alt_intron_to_primary(
                variant, validator
            )

        assert result is False
        assert variant.primary_assembly == "GRCh38"
        assert variant.hgvs_genomic is expected_genomic


    def test_no_current_mapping_uses_other_assembly(
            self, alt_intron_fixture
    ):
        variant, validator = alt_intron_fixture

        validator.hdp.get_tx_mapping_options.return_value = [
            [
                "NM_000001.1",
                "NC_000001.11",
                "splign",
            ],
        ]

        validator.hdp.get_tx_exons.side_effect = [
            make_exons("NW_000001.1"),
            make_exons("NC_000001.11"),
        ]

        expected_genomic = make_hgvs_genomic(
            "NC_000001.11"
        )

        validator.vm.t_to_g.return_value = expected_genomic

        with patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.either_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.start_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "seq_data.supported_for_mapping",
            side_effect=lambda ac, assembly:
                ac == "NC_000001.11"
                and assembly == "GRCh38",
        ):
            result = map_alt_intron_to_primary(
                variant, validator
            )

        assert result is False
        assert variant.primary_assembly == "GRCh38"
        assert variant.hgvs_genomic is expected_genomic


    def test_no_primary_intron_positive_orientation(
            self, alt_intron_fixture
    ):
        variant, validator = alt_intron_fixture

        validator.hdp.get_tx_mapping_options.return_value = []

        validator.hdp.get_tx_exons.return_value = (
            make_exons(
                "NW_000001.1",
                strand=1,
            )
        )

        mapped_genomic = make_hgvs_genomic()
        normalized_genomic = make_hgvs_genomic()
        remapped_transcript = make_remapped_transcript()

        validator.vm.t_to_g.return_value = mapped_genomic
        variant.hn.normalize.return_value = normalized_genomic
        validator.vm.g_to_t.return_value = remapped_transcript

        validator.hdp.get_tx_identity_info.return_value = {
            "hgnc": "TEST"
        }
        validator.hdp.get_pro_ac_for_tx_ac.return_value = (
            "NP_000001.1"
        )

        with patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.either_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.start_position_is_intronic",
            return_value=True,
        ):
            with pytest.raises(
                AltPrimaryMappingError,
                match="AltPrimaryMappingError",
            ):
                map_alt_intron_to_primary(
                    variant, validator
                )

        assert variant.output_type_flag == "gene"
        assert variant.gene_symbol == "TEST"
        assert (
            variant.genome_context_intronic_sequence
            is not None
        )

        variant.hn.normalize.assert_called_once_with(
            mapped_genomic
        )
        variant.reverse_normalizer.normalize.assert_not_called()

        validator.vm.g_to_t.assert_called_once_with(
            normalized_genomic,
            "NM_000001.1",
        )

        # unset_hgvs_obj_ref() should have removed the
        # explicit reference from the remapped transcript.
        assert variant.coding is remapped_transcript
        assert variant.coding.posedit.edit.ref == "A"


    def test_no_primary_intron_negative_orientation(
            self, alt_intron_fixture
    ):
        variant, validator = alt_intron_fixture

        validator.hdp.get_tx_mapping_options.return_value = []

        validator.hdp.get_tx_exons.return_value = (
            make_exons(
                "NW_000001.1",
                strand=-1,
            )
        )

        mapped_genomic = make_hgvs_genomic()
        normalized_genomic = make_hgvs_genomic()
        remapped_transcript = make_remapped_transcript()

        validator.vm.t_to_g.return_value = mapped_genomic

        variant.reverse_normalizer.normalize.return_value = (
            normalized_genomic
        )

        validator.vm.g_to_t.return_value = remapped_transcript

        validator.hdp.get_tx_identity_info.return_value = {
            "hgnc": "TEST"
        }
        validator.hdp.get_pro_ac_for_tx_ac.return_value = (
            "NP_000001.1"
        )

        with patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.either_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.start_position_is_intronic",
            return_value=True,
        ):
            with pytest.raises(
                AltPrimaryMappingError,
                match="AltPrimaryMappingError",
            ):
                map_alt_intron_to_primary(
                    variant, validator
                )

        assert variant.output_type_flag == "gene"

        variant.reverse_normalizer.normalize.assert_called_once_with(
            mapped_genomic
        )
        variant.hn.normalize.assert_not_called()

        validator.vm.g_to_t.assert_called_once_with(
            normalized_genomic,
            "NM_000001.1",
        )

        assert variant.coding is remapped_transcript
        assert variant.coding.posedit.edit.ref == "A"


    def test_ng_no_primary_intron_sets_refseqgene_context(
            self, alt_intron_fixture
    ):
        variant, validator = alt_intron_fixture

        variant.genomic_context_ac = "NG_000001.1"

        validator.hdp.get_tx_mapping_options.return_value = []

        validator.hdp.get_tx_exons.return_value = (
            make_exons(
                "NG_000001.1",
                strand=1,
            )
        )

        mapped_genomic = make_hgvs_genomic(
            "NG_000001.1"
        )
        normalized_genomic = make_hgvs_genomic(
            "NG_000001.1"
        )
        remapped_transcript = make_remapped_transcript()

        validator.vm.t_to_g.return_value = mapped_genomic
        variant.hn.normalize.return_value = normalized_genomic
        validator.vm.g_to_t.return_value = remapped_transcript

        validator.hdp.get_tx_identity_info.return_value = {
            "hgnc": "TEST"
        }
        validator.hdp.get_pro_ac_for_tx_ac.return_value = (
            "NP_000001.1"
        )

        with patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.either_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.start_position_is_intronic",
            return_value=True,
        ):
            with pytest.raises(
                AltPrimaryMappingError
            ):
                map_alt_intron_to_primary(
                    variant, validator
                )

        assert variant.output_type_flag == "gene"

        assert (
            variant.refseqgene_context_intronic_sequence
            is not None
        )

        assert (
            variant.genome_context_intronic_sequence
            is None
        )

        assert variant.coding is remapped_transcript
        assert variant.coding.posedit.edit.ref == "A"


    def test_negative_intronic_boundary(
            self, alt_intron_fixture
    ):
        variant, validator = alt_intron_fixture

        variant.quibble.posedit.pos.start.base = 101
        variant.quibble.posedit.pos.start.offset = -1

        validator.hdp.get_tx_mapping_options.return_value = [
            [
                "NM_000001.1",
                "NC_000001.10",
                "splign",
            ],
        ]

        validator.hdp.get_tx_exons.side_effect = [
            make_exons("NW_000001.1"),
            make_exons("NC_000001.10"),
        ]

        expected_genomic = make_hgvs_genomic(
            "NC_000001.10"
        )

        validator.vm.t_to_g.return_value = expected_genomic

        with patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.either_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.start_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "seq_data.supported_for_mapping",
            side_effect=lambda ac, assembly:
                ac == "NC_000001.10"
                and assembly == "GRCh37",
        ):
            result = map_alt_intron_to_primary(
                variant, validator
            )

        assert result is False
        assert variant.hgvs_genomic is expected_genomic


    def test_intronic_end_position(
            self, alt_intron_fixture
    ):
        variant, validator = alt_intron_fixture

        variant.quibble.posedit.pos.end.base = 100
        variant.quibble.posedit.pos.end.offset = 1

        validator.hdp.get_tx_mapping_options.return_value = [
            [
                "NM_000001.1",
                "NC_000001.10",
                "splign",
            ],
        ]

        validator.hdp.get_tx_exons.side_effect = [
            make_exons("NW_000001.1"),
            make_exons("NC_000001.10"),
        ]

        expected_genomic = make_hgvs_genomic(
            "NC_000001.10"
        )

        validator.vm.t_to_g.return_value = expected_genomic

        with patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.either_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.start_position_is_intronic",
            return_value=False,
        ), patch(
            "VariantValidator.modules.format_converters."
            "seq_data.supported_for_mapping",
            side_effect=lambda ac, assembly:
                ac == "NC_000001.10"
                and assembly == "GRCh37",
        ):
            result = map_alt_intron_to_primary(
                variant, validator
            )

        assert result is False


    def test_identity_sets_genomic_alt_to_ref(
            self, alt_intron_fixture
    ):
        variant, validator = alt_intron_fixture

        variant.quibble.posedit.edit.type = "identity"

        validator.hdp.get_tx_mapping_options.return_value = []
        validator.hdp.get_tx_exons.return_value = (
            make_exons(
                "NW_000001.1",
                strand=1,
            )
        )

        mapped_genomic = make_hgvs_genomic(
            ref="AC",
            alt=None,
        )

        normalized_genomic = make_hgvs_genomic(
            ref="AC",
            alt="AC",
        )

        remapped_transcript = make_remapped_transcript()

        validator.vm.t_to_g.return_value = mapped_genomic
        variant.hn.normalize.return_value = normalized_genomic
        validator.vm.g_to_t.return_value = remapped_transcript

        validator.hdp.get_tx_identity_info.return_value = {
            "hgnc": "TEST"
        }
        validator.hdp.get_pro_ac_for_tx_ac.return_value = None

        with patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.either_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.start_position_is_intronic",
            return_value=True,
        ):
            with pytest.raises(
                AltPrimaryMappingError
            ):
                map_alt_intron_to_primary(
                    variant, validator
                )

        assert mapped_genomic.posedit.edit.alt == "AC"


    def test_genomic_reference_contains_n(
            self, alt_intron_fixture
    ):
        """
        A valid alternate intron that cannot be represented on a
        primary chromosome must stop if the alternate genomic
        reference itself contains undefined N bases.
        """
        variant, validator = alt_intron_fixture

        validator.hdp.get_tx_mapping_options.return_value = []

        validator.hdp.get_tx_exons.return_value = (
            make_exons(
                "NW_000001.1",
                strand=1,
            )
        )

        mapped_genomic = make_hgvs_genomic(
            ref="ANN",
            alt="G",
        )

        normalized_genomic = make_hgvs_genomic(
            ref="ANN",
            alt="G",
        )

        validator.vm.t_to_g.return_value = mapped_genomic
        variant.hn.normalize.return_value = normalized_genomic

        with patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.either_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.start_position_is_intronic",
            return_value=True,
        ):
            with pytest.raises(
                AltPrimaryMappingError,
                match="UndefinedSequenceError",
            ):
                map_alt_intron_to_primary(
                    variant, validator
                )

        assert variant.output_type_flag == "warning"

        assert len(variant.warnings) == 1

        assert (
            "UndefinedSequenceError"
            in variant.warnings[0]
        )

        # We must abort before mapping the N-containing
        # genomic variant back to the transcript.
        validator.vm.g_to_t.assert_not_called()


    def test_missing_identity_and_protein_data(
            self, alt_intron_fixture
    ):
        variant, validator = alt_intron_fixture

        validator.hdp.get_tx_mapping_options.return_value = []
        validator.hdp.get_tx_exons.return_value = (
            make_exons(
                "NW_000001.1",
                strand=1,
            )
        )

        mapped_genomic = make_hgvs_genomic()
        normalized_genomic = make_hgvs_genomic()
        remapped_transcript = make_remapped_transcript()

        validator.vm.t_to_g.return_value = mapped_genomic
        variant.hn.normalize.return_value = normalized_genomic
        validator.vm.g_to_t.return_value = remapped_transcript

        validator.hdp.get_tx_identity_info.side_effect = (
            vvhgvs.exceptions.HGVSDataNotAvailableError(
                "no identity"
            )
        )

        validator.hdp.get_pro_ac_for_tx_ac.side_effect = (
            vvhgvs.exceptions.HGVSDataNotAvailableError(
                "no protein"
            )
        )

        with patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.either_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.start_position_is_intronic",
            return_value=True,
        ):
            with pytest.raises(
                AltPrimaryMappingError,
                match="AltPrimaryMappingError",
            ):
                map_alt_intron_to_primary(
                    variant, validator
                )

        assert variant.output_type_flag == "gene"
        assert variant.coding is remapped_transcript


    def test_mapping_options_ignore_wrong_method_and_non_nc(
            self, alt_intron_fixture
    ):
        variant, validator = alt_intron_fixture

        validator.hdp.get_tx_mapping_options.return_value = [
            [
                "NM_000001.1",
                "NC_000001.10",
                "blat",
            ],
            [
                "NM_000001.1",
                "NW_999999.1",
                "splign",
            ],
            [
                "NM_000001.1",
                "NC_000001.11",
                "splign",
            ],
        ]

        validator.hdp.get_tx_exons.side_effect = [
            make_exons("NW_000001.1"),
            make_exons("NC_000001.11"),
        ]

        expected_genomic = make_hgvs_genomic(
            "NC_000001.11"
        )

        validator.vm.t_to_g.return_value = expected_genomic

        with patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.either_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.start_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "seq_data.supported_for_mapping",
            side_effect=lambda ac, assembly:
                ac == "NC_000001.11"
                and assembly == "GRCh38",
        ):
            result = map_alt_intron_to_primary(
                variant, validator
            )

        assert result is False
        assert variant.primary_assembly == "GRCh38"


    def test_current_primary_mapping_error_is_tolerated(
            self, alt_intron_fixture
    ):
        variant, validator = alt_intron_fixture

        validator.hdp.get_tx_mapping_options.return_value = [
            [
                "NM_000001.1",
                "NC_000001.10",
                "splign",
            ],
        ]

        validator.hdp.get_tx_exons.side_effect = [
            make_exons("NW_000001.1"),
            make_exons("NC_000001.10"),
        ]

        validator.vm.t_to_g.side_effect = (
            vvhgvs.exceptions.HGVSError("mapping failed")
        )

        with patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.either_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.start_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "seq_data.supported_for_mapping",
            side_effect=lambda ac, assembly:
                ac == "NC_000001.10"
                and assembly == "GRCh37",
        ):
            result = map_alt_intron_to_primary(
                variant, validator
            )

        assert result is False


    def test_other_primary_mapping_error_is_tolerated(
            self, alt_intron_fixture
    ):
        variant, validator = alt_intron_fixture

        validator.hdp.get_tx_mapping_options.return_value = [
            [
                "NM_000001.1",
                "NC_000001.11",
                "splign",
            ],
        ]

        validator.hdp.get_tx_exons.side_effect = [
            make_exons("NW_000001.1"),
            make_exons("NC_000001.11"),
        ]

        validator.vm.t_to_g.side_effect = (
            vvhgvs.exceptions.HGVSError("mapping failed")
        )

        with patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.either_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "hgvs_position_utils.start_position_is_intronic",
            return_value=True,
        ), patch(
            "VariantValidator.modules.format_converters."
            "seq_data.supported_for_mapping",
            side_effect=lambda ac, assembly:
                ac == "NC_000001.11"
                and assembly == "GRCh38",
        ):
            result = map_alt_intron_to_primary(
                variant, validator
            )

        assert result is False
        assert variant.primary_assembly == "GRCh38"

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



