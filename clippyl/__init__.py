
#TODO: what documentation goes here?
import unittest

from clippyl.test_build_readid_db import BuildReadidDBTestCase
from clippyl.test_bedgraph_dump import BedgraphDumpTestCase
from clippyl.test_coverage_graphics import CoverageGraphicsTestCase

#TODO: define a unittest test suite to guide test discovery
#load_tests
#make sure that build_readid will be called first during test discovery
#build_ReadidSQLite_dbs(fq_fp_l)

#docs: https://docs.python.org/3/library/unittest.html#load-tests-protocol
test_cases = (BuildReadidDBTestCase, 
              BedgraphDumpTestCase, 
              CoverageGraphicsTestCase)

def load_tests(loader, tests, pattern):
    suite = unittest.TestSuite()
    for test_class in test_cases:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    return suite
