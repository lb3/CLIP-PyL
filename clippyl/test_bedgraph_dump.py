"""test bedgraph-dump on sample data"""
import os

from clippyl.bedgraph_dump import hitsclip_bg_dump
from clippyl.sample_data.paths import (hitsclip_fq_fp_l,
                                       hitsclip_cleav_db_l,
                                       hitsclip_bam_fp_l_discardUnclipped,
                                       hitsclip_bam_fp_l,
                                       hitsclip_n_mapped_reads_l,
                                      )


if __name__ == '__main__':
    
    hitsclip_fq_fp_l()
    hitsclip_cleav_db_l()
    #TODO:check if build_readid will be called first during test discovery
    #build_ReadidSQLite_dbs(fq_fp_l)
    hitsclip_bam_fp_l_discardUnclipped()
    hitsclip_bam_fp_l()
    hitsclip_n_mapped_reads_l()
    
    hitsclip_bg_dump(hitsclip_fq_fp_l(),
                     readid_db_fp_l = hitsclip_cleav_db_l())
#TODO: parameterize the normalization factor because histonly files do not contain all reads

#import random
#import unittest

#class TestSequenceFunctions(unittest.TestCase):

#    def setUp(self):
#        self.seq = list(range(10))

#    def test_shuffle(self):
#        # make sure the shuffled sequence does not lose any elements
#        random.shuffle(self.seq)
#        self.seq.sort()
#        self.assertEqual(self.seq, list(range(10)))

#        # should raise an exception for an immutable sequence
#        self.assertRaises(TypeError, random.shuffle, (1,2,3))

#    def test_choice(self):
#        element = random.choice(self.seq)
#        self.assertTrue(element in self.seq)

#    def test_sample(self):
#        with self.assertRaises(ValueError):
#            random.sample(self.seq, 20)
#        for element in random.sample(self.seq, 5):
#            self.assertTrue(element in self.seq)

#if __name__ == '__main__':
#    unittest.main()

