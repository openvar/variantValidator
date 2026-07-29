from .modules import vvMixinCore as vvMixinCore


class Validator(vvMixinCore.Mixin):
    """
    #Mixins are used to split this very large, complex object over multiple files.
    #There is a logical chain to it, though:
    # vvMixinInit
    #     v
    # vvMixinConverters
    #     v
    # vvMixinCore
    #     v
    # Validator    <- this object.
    """
    pass


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
