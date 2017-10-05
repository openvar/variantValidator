import os
import configuration
import re

# Set App root environment variable
VV_APP_ROOT = os.path.dirname(os.path.abspath(__file__))
os.environ['VV_APP_ROOT'] = VV_APP_ROOT

# Set the version number
fo = open(os.path.join(VV_APP_ROOT, 'VERSION.txt'), 'r')
__version__ = fo.read()
if re.match("^\d+\.\d+\.\d+$", __version__) is not None:
	_is_released_version = True
