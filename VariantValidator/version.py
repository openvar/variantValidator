import pkg_resources
import re
import warnings

# Pull in use_scm_version=True enabled version number
_is_released_version = False
try:
    __version__ = pkg_resources.get_distribution("VariantValidator").version
    if re.match(r"^\d+\.\d+\.\d+$", __version__) is not None:
        _is_released_version = True
except pkg_resources.DistributionNotFound as e:
    warnings.warn("can't get __version__ because %s package isn't installed" % __package__, Warning)
    __version__ = None

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
