import unittest
from unittest.mock import MagicMock, patch

import VariantFormatter.simpleVariantFormatter as svf


class TestSimpleVariantFormatterObjectAPI(unittest.TestCase):

    @patch.object(svf.VariantValidator, "Validator")
    def test_initialises_own_validator(self, mock_validator_class):
        mock_validator = MagicMock()
        mock_validator_class.return_value = mock_validator

        formatter = svf.SimpleVariantFormatter()

        mock_validator_class.assert_called_once_with()
        self.assertIs(formatter.validator, mock_validator)
        self.assertFalse(formatter.testing)

    @patch.object(svf.VariantValidator, "Validator")
    def test_initialises_with_testing_enabled(self, mock_validator_class):
        mock_validator = MagicMock()
        mock_validator_class.return_value = mock_validator

        formatter = svf.SimpleVariantFormatter(testing=True)

        mock_validator_class.assert_called_once_with()
        self.assertIs(formatter.validator, mock_validator)
        self.assertTrue(formatter.testing)

    @patch.object(svf, "_format_impl")
    @patch.object(svf.VariantValidator, "Validator")
    def test_format_passes_all_arguments(
        self,
        mock_validator_class,
        mock_format_impl,
    ):
        mock_validator = MagicMock()
        mock_validator_class.return_value = mock_validator
        mock_format_impl.return_value = {"formatted": "result"}

        formatter = svf.SimpleVariantFormatter(testing=True)

        result = formatter.format(
            variant="NM_000088.4:c.589G>T",
            genome="GRCh38",
            transcript_model="refseq",
            select_transcripts="mane",
            checkOnly=True,
            liftover_level=False,
            legacy_genomic_structure=False,
        )

        self.assertEqual(result, {"formatted": "result"})

        mock_format_impl.assert_called_once_with(
            batch_input="NM_000088.4:c.589G>T",
            genome_build="GRCh38",
            transcript_model="refseq",
            specify_transcripts="mane",
            check_only=True,
            liftover=False,
            validator=mock_validator,
            testing=True,
            legacy_genomic_structure=False,
        )

    @patch.object(svf, "_format_impl")
    @patch.object(svf.VariantValidator, "Validator")
    def test_format_uses_defaults(
        self,
        mock_validator_class,
        mock_format_impl,
    ):
        mock_validator = MagicMock()
        mock_validator_class.return_value = mock_validator
        mock_format_impl.return_value = {}

        formatter = svf.SimpleVariantFormatter()

        formatter.format(
            variant="NM_000088.4:c.589G>T",
            genome="GRCh38",
        )

        mock_format_impl.assert_called_once_with(
            batch_input="NM_000088.4:c.589G>T",
            genome_build="GRCh38",
            transcript_model=None,
            specify_transcripts=None,
            check_only=False,
            liftover=True,
            validator=mock_validator,
            testing=False,
            legacy_genomic_structure=True,
        )

    @patch.object(svf, "_format_impl")
    @patch.object(svf.VariantValidator, "Validator")
    def test_reuses_same_validator_between_calls(
        self,
        mock_validator_class,
        mock_format_impl,
    ):
        mock_validator = MagicMock()
        mock_validator_class.return_value = mock_validator
        mock_format_impl.return_value = {}

        formatter = svf.SimpleVariantFormatter()

        formatter.format(
            variant="NM_000088.4:c.589G>T",
            genome="GRCh38",
        )
        formatter.format(
            variant="NM_000088.4:c.590G>T",
            genome="GRCh38",
        )

        mock_validator_class.assert_called_once_with()
        self.assertEqual(mock_format_impl.call_count, 2)

        for call in mock_format_impl.call_args_list:
            self.assertIs(
                call.kwargs["validator"],
                mock_validator,
            )

    @patch.object(svf, "_format_impl")
    @patch.object(svf.VariantValidator, "Validator")
    def test_string_true_liftover(
        self,
        mock_validator_class,
        mock_format_impl,
    ):
        mock_format_impl.return_value = {}

        formatter = svf.SimpleVariantFormatter()

        formatter.format(
            variant="NM_000088.4:c.589G>T",
            genome="GRCh38",
            liftover_level="True",
        )

        self.assertIs(
            mock_format_impl.call_args.kwargs["liftover"],
            True,
        )

    @patch.object(svf, "_format_impl")
    @patch.object(svf.VariantValidator, "Validator")
    def test_string_false_liftover(
        self,
        mock_validator_class,
        mock_format_impl,
    ):
        mock_format_impl.return_value = {}

        formatter = svf.SimpleVariantFormatter()

        formatter.format(
            variant="NM_000088.4:c.589G>T",
            genome="GRCh38",
            liftover_level="False",
        )

        self.assertIs(
            mock_format_impl.call_args.kwargs["liftover"],
            False,
        )

    @patch.object(svf, "_format_impl")
    @patch.object(svf.VariantValidator, "Validator")
    def test_none_liftover(
        self,
        mock_validator_class,
        mock_format_impl,
    ):
        mock_format_impl.return_value = {}

        formatter = svf.SimpleVariantFormatter()

        formatter.format(
            variant="NM_000088.4:c.589G>T",
            genome="GRCh38",
            liftover_level=None,
        )

        self.assertIs(
            mock_format_impl.call_args.kwargs["liftover"],
            False,
        )

    @patch.object(svf, "_format_impl")
    @patch.object(svf.VariantValidator, "Validator")
    def test_primary_liftover(
        self,
        mock_validator_class,
        mock_format_impl,
    ):
        mock_format_impl.return_value = {}

        formatter = svf.SimpleVariantFormatter()

        formatter.format(
            variant="NM_000088.4:c.589G>T",
            genome="GRCh38",
            liftover_level="primary",
        )

        self.assertEqual(
            mock_format_impl.call_args.kwargs["liftover"],
            "primary",
        )

    @patch.object(svf, "_format_impl")
    @patch.object(svf.VariantValidator, "Validator")
    def test_missing_variant_raises(
        self,
        mock_validator_class,
        mock_format_impl,
    ):
        formatter = svf.SimpleVariantFormatter()

        with self.assertRaisesRegex(
            svf.FormatterSubmissionError,
            "Variant is required",
        ):
            formatter.format(
                variant=None,
                genome="GRCh38",
            )

        mock_format_impl.assert_not_called()

    @patch.object(svf, "_format_impl")
    @patch.object(svf.VariantValidator, "Validator")
    def test_missing_genome_raises(
        self,
        mock_validator_class,
        mock_format_impl,
    ):
        formatter = svf.SimpleVariantFormatter()

        with self.assertRaisesRegex(
            svf.FormatterSubmissionError,
            "Genome build is required",
        ):
            formatter.format(
                variant="NM_000088.4:c.589G>T",
                genome=None,
            )

        mock_format_impl.assert_not_called()

    @patch.object(svf, "_format_impl")
    @patch.object(svf.VariantValidator, "Validator")
    def test_invalid_liftover_raises(
        self,
        mock_validator_class,
        mock_format_impl,
    ):
        formatter = svf.SimpleVariantFormatter()

        with self.assertRaisesRegex(
            svf.FormatterSubmissionError,
            "liftover_level 'invalid' is not supported",
        ):
            formatter.format(
                variant="NM_000088.4:c.589G>T",
                genome="GRCh38",
                liftover_level="invalid",
            )

        mock_format_impl.assert_not_called()


class TestSimpleVariantFormatterFunctionalAPI(unittest.TestCase):

    def setUp(self):
        self.original_global_vfo = svf._GLOBAL_VFO
        self.original_metadata = svf._METADATA
        svf._GLOBAL_VFO = None
        svf._METADATA = None

    def tearDown(self):
        svf._GLOBAL_VFO = self.original_global_vfo
        svf._METADATA = self.original_metadata

    @patch.object(svf.VariantValidator, "Validator")
    def test_global_validator_is_lazy(self, mock_validator_class):
        self.assertIsNone(svf._GLOBAL_VFO)
        mock_validator_class.assert_not_called()

        mock_validator = MagicMock()
        mock_validator_class.return_value = mock_validator

        result = svf._get_global_validator()

        mock_validator_class.assert_called_once_with()
        self.assertIs(result, mock_validator)
        self.assertIs(svf._GLOBAL_VFO, mock_validator)

    @patch.object(svf.VariantValidator, "Validator")
    def test_global_validator_is_reused(self, mock_validator_class):
        mock_validator = MagicMock()
        mock_validator_class.return_value = mock_validator

        first = svf._get_global_validator()
        second = svf._get_global_validator()

        mock_validator_class.assert_called_once_with()
        self.assertIs(first, second)
        self.assertIs(first, mock_validator)

    @patch.object(svf, "_format_impl")
    @patch.object(svf.VariantValidator, "Validator")
    def test_format_creates_validator_lazily(
        self,
        mock_validator_class,
        mock_format_impl,
    ):
        mock_validator = MagicMock()
        mock_validator_class.return_value = mock_validator
        mock_format_impl.return_value = {"formatted": "result"}

        self.assertIsNone(svf._GLOBAL_VFO)

        result = svf.format(
            variant="NM_000088.4:c.589G>T",
            genome="GRCh38",
        )

        self.assertEqual(result, {"formatted": "result"})
        mock_validator_class.assert_called_once_with()
        self.assertIs(svf._GLOBAL_VFO, mock_validator)

        mock_format_impl.assert_called_once_with(
            batch_input="NM_000088.4:c.589G>T",
            genome_build="GRCh38",
            transcript_model=None,
            specify_transcripts=None,
            check_only=False,
            liftover=True,
            validator=mock_validator,
            testing=None,
            legacy_genomic_structure=True,
        )

    @patch.object(svf, "_format_impl")
    @patch.object(svf.VariantValidator, "Validator")
    def test_repeated_format_calls_reuse_global_validator(
        self,
        mock_validator_class,
        mock_format_impl,
    ):
        mock_validator = MagicMock()
        mock_validator_class.return_value = mock_validator
        mock_format_impl.return_value = {}

        svf.format(
            variant="NM_000088.4:c.589G>T",
            genome="GRCh38",
        )
        svf.format(
            variant="NM_000088.4:c.590G>T",
            genome="GRCh38",
        )

        mock_validator_class.assert_called_once_with()
        self.assertEqual(mock_format_impl.call_count, 2)

        for call in mock_format_impl.call_args_list:
            self.assertIs(
                call.kwargs["validator"],
                mock_validator,
            )

    @patch.object(svf, "_format_impl")
    @patch.object(svf.VariantValidator, "Validator")
    def test_supplied_validator_bypasses_global_validator(
        self,
        mock_validator_class,
        mock_format_impl,
    ):
        supplied_validator = MagicMock()
        mock_format_impl.return_value = {}

        svf.format(
            variant="NM_000088.4:c.589G>T",
            genome="GRCh38",
            validator=supplied_validator,
        )

        mock_validator_class.assert_not_called()
        self.assertIsNone(svf._GLOBAL_VFO)

        mock_format_impl.assert_called_once_with(
            batch_input="NM_000088.4:c.589G>T",
            genome_build="GRCh38",
            transcript_model=None,
            specify_transcripts=None,
            check_only=False,
            liftover=True,
            validator=supplied_validator,
            testing=None,
            legacy_genomic_structure=True,
        )

    @patch.object(svf, "_format_impl")
    def test_format_passes_all_arguments(self, mock_format_impl):
        supplied_validator = MagicMock()
        mock_format_impl.return_value = {"formatted": "result"}

        result = svf.format(
            variant="NM_000088.4:c.589G>T",
            genome="GRCh37",
            transcript_model="refseq",
            select_transcripts="mane",
            checkOnly=True,
            liftover_level="primary",
            validator=supplied_validator,
            testing=True,
            legacy_genomic_structure=False,
        )

        self.assertEqual(result, {"formatted": "result"})

        mock_format_impl.assert_called_once_with(
            batch_input="NM_000088.4:c.589G>T",
            genome_build="GRCh37",
            transcript_model="refseq",
            specify_transcripts="mane",
            check_only=True,
            liftover="primary",
            validator=supplied_validator,
            testing=True,
            legacy_genomic_structure=False,
        )

    @patch.object(svf, "_format_impl")
    def test_string_true_liftover(self, mock_format_impl):
        mock_format_impl.return_value = {}
        validator = MagicMock()

        svf.format(
            variant="NM_000088.4:c.589G>T",
            genome="GRCh38",
            liftover_level="True",
            validator=validator,
        )

        self.assertIs(
            mock_format_impl.call_args.kwargs["liftover"],
            True,
        )

    @patch.object(svf, "_format_impl")
    def test_string_false_liftover(self, mock_format_impl):
        mock_format_impl.return_value = {}
        validator = MagicMock()

        svf.format(
            variant="NM_000088.4:c.589G>T",
            genome="GRCh38",
            liftover_level="False",
            validator=validator,
        )

        self.assertIs(
            mock_format_impl.call_args.kwargs["liftover"],
            False,
        )

    @patch.object(svf, "_format_impl")
    def test_none_liftover(self, mock_format_impl):
        mock_format_impl.return_value = {}
        validator = MagicMock()

        svf.format(
            variant="NM_000088.4:c.589G>T",
            genome="GRCh38",
            liftover_level=None,
            validator=validator,
        )

        self.assertIs(
            mock_format_impl.call_args.kwargs["liftover"],
            False,
        )

    @patch.object(svf, "_format_impl")
    def test_missing_variant_raises_before_formatting(
        self,
        mock_format_impl,
    ):
        validator = MagicMock()

        with self.assertRaisesRegex(
            svf.FormatterSubmissionError,
            "Variant is required",
        ):
            svf.format(
                variant=None,
                genome="GRCh38",
                validator=validator,
            )

        mock_format_impl.assert_not_called()

    @patch.object(svf, "_format_impl")
    def test_missing_genome_raises_before_formatting(
        self,
        mock_format_impl,
    ):
        validator = MagicMock()

        with self.assertRaisesRegex(
            svf.FormatterSubmissionError,
            "Genome build is required",
        ):
            svf.format(
                variant="NM_000088.4:c.589G>T",
                genome=None,
                validator=validator,
            )

        mock_format_impl.assert_not_called()

    @patch.object(svf, "_format_impl")
    def test_invalid_liftover_raises_before_formatting(
        self,
        mock_format_impl,
    ):
        validator = MagicMock()

        with self.assertRaisesRegex(
            svf.FormatterSubmissionError,
            "liftover_level 'invalid' is not supported",
        ):
            svf.format(
                variant="NM_000088.4:c.589G>T",
                genome="GRCh38",
                liftover_level="invalid",
                validator=validator,
            )

        mock_format_impl.assert_not_called()


class TestSimpleVariantFormatterSharedHelpers(unittest.TestCase):

    def test_normalise_transcript_selection_all(self):
        self.assertEqual(
            svf._normalise_transcript_selection('["all"]'),
            "all",
        )

    def test_normalise_transcript_selection_raw(self):
        self.assertEqual(
            svf._normalise_transcript_selection('["raw"]'),
            "raw",
        )

    def test_normalise_transcript_selection_mane(self):
        self.assertEqual(
            svf._normalise_transcript_selection('["mane"]'),
            "mane",
        )

    def test_normalise_transcript_selection_mane_select(self):
        self.assertEqual(
            svf._normalise_transcript_selection('["mane_select"]'),
            "mane_select",
        )

    def test_normalise_transcript_selection_select(self):
        self.assertEqual(
            svf._normalise_transcript_selection('["select"]'),
            "select",
        )

    def test_normalise_transcript_selection_preserves_other_values(self):
        self.assertEqual(
            svf._normalise_transcript_selection("NM_000088.4"),
            "NM_000088.4",
        )

    def test_normalise_batch_input_list(self):
        variants = [
            "NM_000088.4:c.589G>T",
            "NM_000088.4:c.590G>T",
        ]

        self.assertEqual(
            svf._normalise_batch_input(variants),
            variants,
        )

    def test_normalise_batch_input_json_list(self):
        result = svf._normalise_batch_input(
            '["NM_000088.4:c.589G>T", "NM_000088.4:c.590G>T"]'
        )

        self.assertEqual(
            result,
            [
                "NM_000088.4:c.589G>T",
                "NM_000088.4:c.590G>T",
            ],
        )

    def test_normalise_batch_input_single_variant(self):
        self.assertEqual(
            svf._normalise_batch_input(
                "NM_000088.4:c.589G>T"
            ),
            ["NM_000088.4:c.589G>T"],
        )

    def test_contains_hgvs_type(self):
        self.assertTrue(
            svf._contains_hgvs_type(
                "NM_000088.4:c.589G>T"
            )
        )

    def test_does_not_contain_hgvs_type(self):
        self.assertFalse(
            svf._contains_hgvs_type(
                "17-50198002-C-A"
            )
        )


if __name__ == "__main__":
    unittest.main()


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
