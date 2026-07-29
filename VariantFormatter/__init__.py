import importlib.metadata
import re


__version__ = importlib.metadata.version("VariantFormatter")
_is_released_version = re.fullmatch(r"\d+\.\d+\.\d+", __version__) is not None


# Copyright (C) 2016-2026 VariantValidator Contributors
# This file is part of VariantValidator and is distributed under the
# GNU Affero General Public License, version 3 or (at your option) any
# later version. See the LICENSE file in the project root for the full
# licence terms.
# SPDX-License-Identifier: AGPL-3.0-or-later
