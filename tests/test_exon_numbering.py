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
        cls.vv.testing = True

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
        results = results['NM_000088.3:c.642=']['variant_exonic_positions']
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

    def test_8(self):
        results = self.vv.validate("NM_004656.3:c.2189_*2del", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_004656.3:c.2189_*2del']['variant_exonic_positions']
        print(results)
        assert results == {
            "NC_000003.11": {
                "end_exon": "17",
                "start_exon": "17"
            },
            "NC_000003.12": {
                "end_exon": "17",
                "start_exon": "17"
            }
        }

    def test_8a(self):
        results = self.vv.validate("NM_004656.3:c.*1_*2del", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_004656.3:c.*1_*2del']['variant_exonic_positions']
        print(results)
        assert results == {
            "NC_000003.11": {
                "end_exon": "17",
                "start_exon": "17"
            },
            "NC_000003.12": {
                "end_exon": "17",
                "start_exon": "17"
            }
        }

    def test_9(self):
        results = self.vv.validate("NM_004006.2:c.-244T>A", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_004006.2:c.-244T>A']['variant_exonic_positions']
        print(results)
        assert results == {
            "NC_000023.10": {
                "end_exon": "1",
                "start_exon": "1"
            },
            "NC_000023.11": {
                "end_exon": "1",
                "start_exon": "1"
            },
            "NG_012232.1": {
                "end_exon": "1",
                "start_exon": "1"
            }
        }

    def test_9a(self):
        results = self.vv.validate("NM_004006.2:c.-244_-200del", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_004006.2:c.-242_-198del']['variant_exonic_positions']
        print(results)
        assert results == {
            "NC_000023.10": {
                "end_exon": "1",
                "start_exon": "1"
            },
            "NC_000023.11": {
                "end_exon": "1",
                "start_exon": "1"
            },
            "NG_012232.1": {
                "end_exon": "1",
                "start_exon": "1"
            }
        }

    def test_9c(self):
        results = self.vv.validate("NM_004006.2:c.-2_3del", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_004006.2:c.-2_3del']['variant_exonic_positions']
        print(results)
        assert results == {
            "NC_000023.10": {
                "end_exon": "1",
                "start_exon": "1"
            },
            "NC_000023.11": {
                "end_exon": "1",
                "start_exon": "1"
            },
            "NG_012232.1": {
                "end_exon": "1",
                "start_exon": "1"
            }
        }

    def test_10(self):
        results = self.vv.validate("NM_001362177.1:c.-1+1dup", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_001362177.1:c.-1+1dup']['variant_exonic_positions']
        print(results)
        assert results == {
            "NC_000009.11": {
                "end_exon": "4i",
                "start_exon": "4i"
            },
            "NC_000009.12": {
                "end_exon": "4i",
                "start_exon": "4i"
            }
        }

    def test_10a(self):
        results = self.vv.validate("NM_001354304.1:c.-95-121A>G", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_001354304.1:c.-95-121A>G']['variant_exonic_positions']
        print(results)
        assert results == {
            "NC_000012.11": {
                "end_exon": "1i",
                "start_exon": "1i"
            },
            "NC_000012.12": {
                "end_exon": "1i",
                "start_exon": "1i"
            }
        }

    def test_10b(self):
        results = self.vv.validate("NM_001362177.1:c.-1+1_1-1del", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_001362177.1:c.-1+1_1-1del']['variant_exonic_positions']
        print(results)
        assert results == {
            "NC_000009.11": {
                "end_exon": "4i",
                "start_exon": "4i"
            },
            "NC_000009.12": {
                "end_exon": "4i",
                "start_exon": "4i"
            }
        }

    def test_10c(self):
        results = self.vv.validate("NM_001354304.1:c.-96+5_-95-5A>G", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_001354304.1:c.-96+6_-95-5del']['variant_exonic_positions']
        print(results)
        assert results == {
            "NC_000012.11": {
                "end_exon": "1i",
                "start_exon": "1i"
            },
            "NC_000012.12": {
                "end_exon": "1i",
                "start_exon": "1i"
            }
        }

    def test_11(self):
        results = self.vv.validate("NM_058197.4:c.*74-1G>T", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_058197.4:c.*74-1G>T']['variant_exonic_positions']
        print(results)
        assert results == {
            "NC_000009.11": {
                "end_exon": "1i",
                "start_exon": "1i"
            },
            "NC_000009.12": {
                "end_exon": "1i",
                "start_exon": "1i"
            }
        }

    def test_11a(self):
        results = self.vv.validate("NM_058197.4:c.*73+12_*74-12del", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_058197.4:c.*73+12_*74-12del']['variant_exonic_positions']
        print(results)
        assert results == {
            "NC_000009.11": {
                "end_exon": "1i",
                "start_exon": "1i"
            },
            "NC_000009.12": {
                "end_exon": "1i",
                "start_exon": "1i"
            }
        }

    def test_11b(self):
        results = self.vv.validate("NM_058197.4:c.*73+12=", 'GRCh38', 'all').format_as_dict(test=True)
        results = results['NM_058197.4:c.*73+12=']['variant_exonic_positions']
        print(results)
        assert results == {
            "NC_000009.11": {
                "end_exon": "1i",
                "start_exon": "1i"
            },
            "NC_000009.12": {
                "end_exon": "1i",
                "start_exon": "1i"
            }
        }


if __name__ == "__main__":
    unittest.main()


# <LICENSE>
# Copyright (C) 2016-2024 VariantValidator Contributors
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
