"""test bedgraph-dump on sample data"""
import os

from clippyl.sample_data.paths import hitsclip_bam_dir
from clippyl.sample_data.paths import hitsclip_discardUnclipped_fq_dir
from clippyl.bedgraph_dump import hitsclip_bg_dump

if __name__ == '__main__':
    # getting relevant filepaths lists
    fp_l = os.listdir(hitsclip_bam_dir())
    fp_l = [os.path.join(hitsclip_bam_dir(), fp) for fp in fp_l]
    bam_fp_l = [fp for fp in fp_l if fp.split('.')[-1] == 'bam']
    bam_fp_l_discardUnclipped = sorted([fp for fp in bam_fp_l if fp.split('.')[2] == 'discardUnclipped'])
    bam_fp_l = sorted([fp for fp in bam_fp_l if fp.split('.')[2] != 'discardUnclipped'])
    
    #NOTE: run test_build_readid_db must to generate readid databases
    fp_l = os.listdir(hitsclip_discardUnclipped_fq_dir())
    fp_l = [os.path.join(hitsclip_discardUnclipped_fq_dir(), fp) for fp in fp_l]
    readid_db_fp_l = sorted([fp for fp in fp_l if fp.split('.')[-1] == 'readids'])
    
    print (readid_db_fp_l)
    
    #hitsclip_bg_dump(bam_fp_l_discardUnclipped,)
    
    hitsclip_bg_dump(bam_fp_l, readid_db_fp_l = readid_db_fp_l)
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

