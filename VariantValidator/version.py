import importlib.metadata
import re
import warnings

# Pull in use_scm_version=True enabled version number
_is_released_version = False
try:
    __version__ = importlib.metadata.version("VariantValidator")
    if re.match(r"^\d+\.\d+\.\d+$", __version__) is not None:
        _is_released_version = True
except importlib.metadata.PackageNotFoundError:
    warnings.warn("can't get __version__ because VariantValidator package isn't installed", Warning)
    __version__ = None

# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
