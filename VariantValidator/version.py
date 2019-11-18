import pkg_resources
import re
import warnings

# Pull in use_scm_version=True enabled version number
_is_released_version = False
try:
    __version__ = pkg_resources.get_distribution("rest_VariantValidator").version
    if re.match(r"^\d+\.\d+\.\d+$", __version__) is not None:
        _is_released_version = True
except pkg_resources.DistributionNotFound as e:
    warnings.warn("can't get __version__ because %s package isn't installed" % __package__, Warning)
    __version__ = None
