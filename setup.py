#!/usr/bin/env python3

import sys
import os.path
from distutils.core import setup


if sys.version_info[0] >= 3 and sys.version_info[1] > 2:
   sys.stderr.write( "Error in setup script for CLIP-PyL:\n" )
   sys.stderr.write( "You need at least version 3.2 of Python to use CLIP-PyL.\n" )
   sys.exit( 1 )

#TODO: check for pysam and matplotlib dependencies with try: import, except

# docs:
# https://docs.python.org/3/distutils/index.html#distutils-index
setup(
    name='CLIP-PyL',
    version='0.0.1-alpha',
    author='Lionel Brooks 3rd, PhD',
    author_email='lbrooks.thethird+cp@gmail.com',
    url='https://github.com/lb3/CLIP-PyL',
    packages=['clippyl'],
    package_data={'clippyl': ['sample_data/*',
                              'sample_data/cis_element_annotations/*',
                              'sample_data/genomic_interval_queries/*',
                              'sample_data/HITS-CLIP_SLBP_histone_mRNA_01/adapter_clipped_reads_fastq/*',
                              'sample_data/HITS-CLIP_SLBP_histone_mRNA_01/bwa_samse_hg19/*'],
                 },
    scripts=['scripts/clippyl-bedgraph',
             'scripts/clippyl-build-readid-db',
             'scripts/clippyl-graphics'],
    license='MIT License', #SEE LICENSE.txt file included with this package
    description='tools for finding crosslinked nucleotides in aligned CLIP-seq data',
    long_description=open('README').read()
)
