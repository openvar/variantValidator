#!/usr/bin/env python

# Prefer setuptools over distutils
from setuptools import setup, find_packages

setup(
    name='VariantValidator',
    version='0.1.0',
    description='API for accurate, mapping and formatting of sequence variants using HGVS nomenclature',
    long_description=open('README.txt').read(),
    url='',
    author='Peter J. Causey-Freeman',
    author_email='pjf9@leicester.ac.uk',
    license=open('LICENSE.txt').read(),
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

        # License
        license=open('LICENSE.txt').read(),

        # Specify the Python versions
        'Programming Language :: Python :: 2.7',
    ],
 
    # What does your project relate to?
    keywords='sample setuptools development',

    packages=['variantValidator', ],
	# List run-time dependencies here.  These will be installed by pip when the project is installed.
    install_requires=[
        "hgvs >= 1.0.0", # This will install BioPython
		"biocommons.seqrepo >= 0.3.5",
		"httplib2 >= 0.9.0",
		"configparser >= 3.5.0",
		"pyliftover >= 0.3",
		"biotools >= 0.3.0",
		"mysql_connector >=	2.1.4",  
    ],
)