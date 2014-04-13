#!/usr/bin/env python

import sys
import os.path

try:
   from setuptools import setup, Extension
except ImportError:   
   sys.stderr.write( "Could not import 'setuptools', falling back to 'distutils'.\n" )
   from distutils.core import setup, Extension

if sys.version_info[0] < 2 or sys.version_info < 5:
   sys.stderr.write( "Error in setup script for HTSeq:\n" )
   sys.stderr.write( "You need at least version 2.5 of Python to use HTSeq.\n" )
   sys.exit( 1 )

if sys.version_info[0] >= 3:
   sys.stderr.write( "Error in setup script for HTSeq:\n" )
   sys.stderr.write( "Sorry, this package does not yet work with Python 3.\n" )
   sys.stderr.write( "Please use Python 2.x, x>=5.\n" )
   sys.exit( 1 )

try:
   import numpy
except ImportError:
   sys.stderr.write( "Setup script for HTSeq: Failed to import 'numpy'.\n" )
   sys.stderr.write( "Please install numpy and then try again to install HTSeq.\n" )
   sys.exit( 1 )
   
numpy_include_dir = os.path.join( os.path.dirname( numpy.__file__ ), 'core', 'include' )


setup(
    name='CLIP-PyL',
    version='0.1dev',
    packages=['clippyl',],
    scripts=['bin/...'],
    license='LICENSE.txt',
    description='tools for finding crosslinked nucleotides in aligned CLIP-seq data',
    long_description=open('README.txt').read()
)
