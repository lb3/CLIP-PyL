
#TODO: what documentation goes here?
import unittest

from clippyl.test_bedgraph_dump import BG_Dump_Basic_TC
from clippyl.test_coverage_graphics import Coverage_Graphics_Basic_TC

# only the most basic functionality of clippyl is tested here
# run the individual test_* modules for more comprehensive
# tests that include the readid db and cis element sl3 feature

#docs: https://docs.python.org/3/library/unittest.html#load-tests-protocol
test_cases = (BG_Dump_Basic_TC, 
              Coverage_Graphics_Basic_TC)

def load_tests(loader, tests, pattern):
    suite = unittest.TestSuite()
    for test_class in test_cases:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    return suite
