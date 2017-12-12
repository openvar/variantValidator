# -*- coding: utf-8 -*-

"""
Add the location of the config.ini as an environment varaible
"""

import os
import re
CONF_ROOT = os.path.dirname(os.path.abspath(__file__))
os.environ['CONF_ROOT'] = CONF_ROOT


# <LICENSE>

# </LICENSE>
