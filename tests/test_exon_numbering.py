"""
Exon_numbering_tests
Authors: Katie Williams (@kwi11iams) and Katherine Winfield (@kjwinfield)
This code runs tests on the module exon_numbering.py to check the outputs are as expected
"""
import unittest
from VariantValidator import Validator


class TestExonNumbering(unittest.TestCase):

    @classmethod
    def setup_class(cls):
        cls.vv = Validator()

    """
    Class TextExonNumbering automates running the tests, and reports failure if
    the output is not what is expected
    """
    def test_1(self):
        results = self.vv.validate("NM_007294.3:c.1067A>G", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_007294.3:c.1067A>G']['variant_exonic_positions']
        print(results)
        assert results == {
            "NC_000017.10": {
                "end_exon": "10",
                "start_exon": "10"
            },
            "NC_000017.11": {
                "end_exon": "10",
                "start_exon": "10"
            },
            "NG_005905.2": {
                "end_exon": "10",
                "start_exon": "10"
            },
        }

    def test_2(self):
        results = self.vv.validate("NM_000088.3:c.642+1G>A", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_000088.3:c.642+1G>A']['variant_exonic_positions']
        print(results)
        assert results == {
                "NC_000017.10": {
                    "end_exon": "8i",
                    "start_exon": "8i"
                },
                "NC_000017.11": {
                    "end_exon": "8i",
                    "start_exon": "8i"
                },
                "NG_007400.1": {
                    "end_exon": "8i",
                    "start_exon": "8i"
                }
            }

    def test_2b(self):
        results = self.vv.validate("NM_000088.3:c.642=", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_000088.3:c.642A=']['variant_exonic_positions']
        print(results)
        assert results == {
                "NC_000017.10": {
                    "end_exon": "8",
                    "start_exon": "8"
                },
                "NC_000017.11": {
                    "end_exon": "8",
                    "start_exon": "8"
                },
                "NG_007400.1": {
                    "end_exon": "8",
                    "start_exon": "8"
                }
            }

    def test_3(self):
        results = self.vv.validate("NM_000094.3:c.6751-3_6751-2del", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_000094.3:c.6751-3_6751-2del']['variant_exonic_positions']
        print(results)
        assert results == {
                "NC_000003.11": {
                    "end_exon": "85i",
                    "start_exon": "85i"
                },
                "NC_000003.12": {
                    "end_exon": "85i",
                    "start_exon": "85i"
                },
                "NG_007065.1": {
                    "end_exon": "85i",
                    "start_exon": "85i"
                }
            }

    def test_4(self):
        results = self.vv.validate("NM_000088.3:c.589G>T", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_000088.3:c.589G>T']['variant_exonic_positions']
        print(results)
        assert results == {
                "NC_000017.10": {
                    "end_exon": "8",
                    "start_exon": "8"
                },
                "NC_000017.11": {
                    "end_exon": "8",
                    "start_exon": "8"
                },
                "NG_007400.1": {
                    "end_exon": "8",
                    "start_exon": "8"
                }
            }

    def test_4a(self):
        results = self.vv.validate("NM_000088.3:c.589-1G>T", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_000088.3:c.589-1G>T']['variant_exonic_positions']
        print(results)
        assert results == {
                "NC_000017.10": {
                    "end_exon": "7i",
                    "start_exon": "7i"
                },
                "NC_000017.11": {
                    "end_exon": "7i",
                    "start_exon": "7i"
                },
                "NG_007400.1": {
                    "end_exon": "7i",
                    "start_exon": "7i"
                }
            }

    def test_5(self):
        results = self.vv.validate("NM_000088.3:c.642del", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_000088.3:c.642del']['variant_exonic_positions']
        print(results)
        assert results == {
                "NC_000017.10": {
                    "end_exon": "8",
                    "start_exon": "8"
                },
                "NC_000017.11": {
                    "end_exon": "8",
                    "start_exon": "8"
                },
                "NG_007400.1": {
                    "end_exon": "8",
                    "start_exon": "8"
                }
            }

    def test_6(self):
        results = self.vv.validate("NM_000088.3:c.642+5_643-25del", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_000088.3:c.642+5_643-25del']['variant_exonic_positions']
        print(results)
        assert results == {
                "NC_000017.10": {
                    "end_exon": "8i",
                    "start_exon": "8i"
                },
                "NC_000017.11": {
                    "end_exon": "8i",
                    "start_exon": "8i"
                },
                "NG_007400.1": {
                    "end_exon": "8i",
                    "start_exon": "8i"
                }
            }

    def test_6a(self):
        results = self.vv.validate("NM_000088.3:c.642+5_642+25del", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_000088.3:c.642+5_642+25del']['variant_exonic_positions']
        print(results)
        assert results == {
                "NC_000017.10": {
                    "end_exon": "8i",
                    "start_exon": "8i"
                },
                "NC_000017.11": {
                    "end_exon": "8i",
                    "start_exon": "8i"
                },
                "NG_007400.1": {
                    "end_exon": "8i",
                    "start_exon": "8i"
                }
            }

    def test_6b(self):
        results = self.vv.validate("NM_000088.3:c.643-25_643-5del", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_000088.3:c.643-25_643-5del']['variant_exonic_positions']
        print(results)
        assert results == {
                "NC_000017.10": {
                    "end_exon": "8i",
                    "start_exon": "8i"
                },
                "NC_000017.11": {
                    "end_exon": "8i",
                    "start_exon": "8i"
                },
                "NG_007400.1": {
                    "end_exon": "8i",
                    "start_exon": "8i"
                }
            }

    def test_7(self):
        results = self.vv.validate("NM_000094.3:c.6751-3_6753del", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_000094.3:c.6751-3_6753del']['variant_exonic_positions']
        print(results)
        assert results == {
                "NC_000003.11": {
                    "end_exon": "86",
                    "start_exon": "85i"
                },
                "NC_000003.12": {
                    "end_exon": "86",
                    "start_exon": "85i"
                },
                "NG_007065.1": {
                    "end_exon": "86",
                    "start_exon": "85i"
                }
            }


if __name__ == "__main__":
    unittest.main()


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