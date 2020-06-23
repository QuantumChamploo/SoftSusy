#! /usr/bin/env python

from distutils.core import setup

## Setup definition
import mypyslha
__doc__ = pyslha.__doc__

setup(name = 'mypyslha',
      version = mypyslha.__version__,
      py_modules = ["mypyslha"],
      scripts = ["slhaplot", "slha2isawig", "isawig2slha"],
      install_requires = ["tex2pix >=0.1.5"],
      author = 'Andy Buckley',
      author_email = 'andy@insectnation.org',
      url = 'http://www.insectnation.org/projects/pyslha',
      description = 'Parsing, manipulating, and visualising SUSY Les Houches Accord data',
      long_description = __doc__,
      keywords = 'supersymmetry susy slha simulation mass decay hep physics particle',
      license = 'GPL')
