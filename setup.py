from distutils.core import setup

setup(
    name='CLIP-PyL',
    version='0.1dev',
    packages=['clippyl',],
    scripts=['bin/...'],
    license='LICENSE.txt',
    description='tools for finding crosslinked nucleotides in aligned CLIP-seq data',
    long_description=open('README.txt').read()
)
