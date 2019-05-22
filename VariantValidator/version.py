import os

with open(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'VERSION.txt')) as fh:
    __version__ = fh.read().strip()
