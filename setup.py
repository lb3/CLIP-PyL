#!/usr/bin/env python3

import sys
import os.path
from distutils.core import setup


if sys.version_info[0] >= 3 and sys.version_info[1] > 2:
   sys.stderr.write( "Error in setup script for CLIP-PyL:\n" )
   sys.stderr.write( "You need at least version 3 of Python to use CLIP-PyL.\n" )
   sys.exit( 1 )

setup(
    name='CLIP-PyL',
    version='00.01',
    author='Lionel Brooks 3rd, PhD',
    author_email='',
    packages=['clippyl'],
    package_data={'clippyl': ['clippyl/sample_data/*']},
    scripts=['scripts/clippyl-bedgraph',
             'scripts/clippyl-build-readid-db',
             'scripts/clippyl-graphics'],
    license='LICENSE.txt',
    description='tools for finding crosslinked nucleotides in aligned CLIP-seq data',
    long_description=open('README.md').read()
)
