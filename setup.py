#!/usr/bin/env python3

import sys
import os.path
from distutils.core import setup


if sys.version_info[0] != 3 or sys.version_info[1] < 2:
   sys.stderr.write( "Error in setup script for HTSeq:\n" )
   sys.stderr.write( "You need at least version 2.5 of Python to use HTSeq.\n" )
   sys.exit( 1 )

setup(
    name='CLIP-PyL',
    version='0.1dev',
    packages=['clippyl'],
    package_dir={'mypkg': 'clippyl'},
    package_data={'mypkg': ['clippyl/sample_data/*']},
    scripts=['scripts/build_hitsclip_PyL'],
    license='LICENSE.txt',
    description='tools for finding crosslinked nucleotides in aligned CLIP-seq data',
    long_description=open('README.txt').read()
)
