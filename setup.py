#!/usr/bin/env python

# Prefer setuptools over distutils
from setuptools import setup, find_packages

# with open('VariantValidator/version.py') as ins:
#     version = ins.read()
#     version = version.split('=')[1].strip()
#     version = version.replace("'", "")

setup(
    name='VariantValidator',
    description='API for accurate, mapping and formatting of sequence variants using HGVS nomenclature',
    long_description=open('README.md').read(),
    url='https://github.com/openvar/variantValidator',
    use_scm_version=True,
    zip_safe=True,
    author="VariantValidator Contributors",
    author_email='admin@variantvalidator.org',
    packages=['VariantValidator', 'VariantValidator.modules'],
    include_package_data=True,
    license="GNU AFFERO GENERAL PUBLIC LICENSE, Version 3 (https://www.gnu.org/licenses/agpl-3.0.en.html)",
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Audience
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        # Specify the Python versions
        'Programming Language :: Python :: 3.6',
    ],
    scripts=[
        'bin/update_vdb.py',
        'bin/variant_validator.py',
        'bin/vv_configure.py'
    ],
    data_files=[
        ('configuration', ['configuration/default.ini', 'configuration/empty_vv_db.sql'])
    ],
    # What does your project relate to?
    keywords=[
        "bioinformatics",
        "computational biology",
        "genome variants",
        "genome variation",
        "genomic variants",
        "genomic variation",
        "genomics",
        "hgvs",
        "HGVS",
        "sequencevariants",
    ],

    # List run-time dependencies here.  These will be installed by pip when the project is installed.
    install_requires=[
        "biocommons.seqrepo>=0.5.1",
        "httplib2>=0.9.0",
        "configparser>=3.5.0",
        "pyliftover>=0.3",
        "biotools>=0.3.0",
        # This version has been tested for use with vvhgvs (ignore warnings) and needs to be maintained to ensure
        # biopython is still installed when vvhgvs is updated in-line with biocommons.hgvs which has had the dependancy
        # removed
        "biopython==1.74",
        "requests",
        "mysql-connector-python",
    ],
    setup_requires=[
        "setuptools_scm",
    ]

)

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
