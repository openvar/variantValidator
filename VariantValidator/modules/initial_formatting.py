import re
import logging

logger = logging.getLogger(__name__)


def remove_gene_symbol_from_ref(my_variant, validator):
    """
    Remove a redundant gene symbol from the reference accession while
    preserving the HGVS SequenceVariant object.
    """
    accession = my_variant.quibble.ac

    start = accession.find('(')
    if start == -1:
        return

    end = accession.find(')', start + 1)
    if end == -1:
        return

    gene_symbol = accession[start + 1:end]
    if not gene_symbol:
        return

    hgnc_symbol = validator.db.get_hgnc_symbol(gene_symbol)

    if hgnc_symbol == 'none' and not gene_symbol.startswith('MT-'):
        return

    warning = (
        f'VariantSyntaxError: Removing redundant gene symbol '
        f'{hgnc_symbol} from variant description'
    )

    my_variant.quibble.ac = accession[:start] + accession[end + 1:]
    my_variant.warnings.append(warning)
    logger.info(warning)


def initial_user_formattng(my_variant, validator):
    """
    Compile HGVS input errors and provide improved warnings/error messages
    while retaining the parsed HGVS SequenceVariant.
    """
    # Remove redundant gene symbols from the reference accession.
    remove_gene_symbol_from_ref(my_variant, validator)

    # Inspect the original submitted description for redundant reference
    # bases which may no longer be represented after HGVS parsing.
    redundant_bases = re.search(
        r'(del|inv|dup|ins)[GATC]+',
        my_variant.original
    )

    if redundant_bases and redundant_bases.group(1) != 'ins':
        warning = (
            'VariantSyntaxError: Removing redundant reference bases '
            'from variant description'
        )
        my_variant.warnings.append(warning)
        logger.info(warning)

    # All reference sequence identifiers must use the correct case.
    # LRG transcript identifiers are the exception because the transcript
    # separator uses a lowercase t, e.g. LRG_1t1.
    original_ac, _sep, _remainder = my_variant.original.partition(':')
    uc_original_ac = original_ac.upper()

    if uc_original_ac.startswith('LRG'):
        if not original_ac.startswith('LRG'):
            error = (
                'This not a valid HGVS description, due to characters being '
                'in the wrong case. Please check the use of upper- and '
                'lowercase characters.'
            )
            my_variant.warnings.append(error)
            logger.info(error)

        if 'T' in original_ac:
            error = (
                'This not a valid HGVS description, due to characters being '
                'in the wrong case. Please check the use of upper- and '
                'lowercase characters.'
            )
            my_variant.warnings.append(error)
            logger.info(error)

            my_variant.quibble.ac = my_variant.quibble.ac.replace('T', 't')

    elif original_ac != uc_original_ac:
        if uc_original_ac.startswith(('CHR', 'GRCH', 'HG')):
            error = (
                'This is not a valid HGVS variant description, because no '
                'reference sequence ID has been provided'
            )
        else:
            error = (
                'This not a valid HGVS description, due to characters being '
                'in the wrong case. Please check the use of upper- and '
                'lowercase characters.'
            )

        my_variant.warnings.append(error)
        logger.info(error)

    elif uc_original_ac.startswith(('CHR', 'GRCH', 'HG')):
        error = (
            'This is not a valid HGVS variant description, because no '
            'reference sequence ID has been provided'
        )
        my_variant.warnings.append(error)
        logger.info(error)

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
