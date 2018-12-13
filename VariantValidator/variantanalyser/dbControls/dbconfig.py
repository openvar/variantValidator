# -*- coding: utf-8 -*-

"""
dbconfig.py

Configures the MySQL connection using the config.ini or environment variables
"""

from configparser import ConfigParser
import os

# Get the conf_root from the os
CONF_ROOT = os.environ.get('CONF_ROOT')
def read_db_config(filename=os.path.join(CONF_ROOT, 'config.ini'), section='mysql'):
    """ Read database configuration file and return a dictionary object
    :param filename: name of the configuration file
    :param section: section of database configuration
    :return: a dictionary of database parameters
    """
    # create parser and read ini configuration file
    parser = ConfigParser()
    with open(filename) as f:
        parser.read_file(f)
    
    # get section, default to mysql
    db = {}
    if parser.has_section(section):
    	items = parser.items(section)
    	for item in items:
    		db[item[0]] = item[1]
    else:
    	raise Exception('{0} not found in the {1} file'.format(section, filename))
    
    return db
    
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