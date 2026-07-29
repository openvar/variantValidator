
def methyl_syntax(my_variant):
    """
    Remove and store a methylation syntax suffix before HGVS object parsing.

    :param my_variant:
    :return: Updated variant if methylation syntax is detected, otherwise None.
    """
    quibble = my_variant.quibble

    if "|" not in quibble:
        return None

    if "|gom" in quibble:
        my_variant.reformat_output = "|gom"
    elif "|lom" in quibble:
        my_variant.reformat_output = "|lom"
    elif "|met=" in quibble:
        my_variant.reformat_output = "|met="
    else:
        return None

    my_variant.quibble = quibble.split("|", 1)[0] + "="
    return my_variant


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
