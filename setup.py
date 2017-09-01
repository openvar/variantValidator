#!/usr/bin/env python

from distutils.core import setup

setup(
    name='VariantValidator',
    version='0.1.0',
    author='Peter J. Causey-Freeman',
    author_email='pjf9@leicester.ac.uk',
    packages=['variantValidator', ],
    url='',
    license='LICENSE.txt',
    description='Tool for accurate, mapping and formatting of sequence variants using HGVS nomenclature',
    long_description=open('README.txt').read(),
    install_requires=[
        "hgvs >= 1.0.0", # This will install BioPython
		"biocommons.seqrepo >= 0.3.5",
		"httplib2 >= 0.9.0",
		"configparser >= 3.5.0",
		"pyliftover >= 0.3",
		"biotools >= 0.3.0",
		"mysql_connector >=	2.1.4"  
    ],
)