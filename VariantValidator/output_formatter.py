# -*- coding: utf-8 -*-

"""
Module that will be embedded into the valstr method and translate method

Two functions,

format nucleotide descriptions removes the ref bases from displayed descriptions

format protein descriptions return the description using the single letter aa alphabet
"""

import hgvs

"""
format protein description into single letter aa code
"""


def single_letter_protein(hgvs_protein):
    hgvs_protein_slc = hgvs_protein.format({'p_3_letter': False})
    return hgvs_protein_slc


"""
format nucleotide descriptions to not display reference base
"""


def remove_reference(hgvs_nucleotide):
    hgvs_nucleotide_refless = hgvs_nucleotide.format({'max_ref_length': 0})
    return hgvs_nucleotide_refless

# <LICENSE>
# Copyright (C) 2018  Peter Causey-Freeman, University of Leicester
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