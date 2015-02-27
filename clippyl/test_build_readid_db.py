import os
import unittest

from clippyl.build_readid_db import build_ReadidSQLite_dbs
from clippyl.sample_data.paths import hitsclip_discardUnclipped_fq_dir

class BuildReadidDBTestCase(unittest.TestCase):
    
    def runTest(self):
        # getting relevant sample data filepaths
        fp_l = os.listdir(hitsclip_discardUnclipped_fq_dir())
        fp_l = [os.path.join(hitsclip_discardUnclipped_fq_dir(), fp) for fp in fp_l]
        fq_fp_l = [fp for fp in fp_l if fp.split('.')[-2:] == ['fq','gz']]
        print (fq_fp_l)
        
        build_ReadidSQLite_dbs(fq_fp_l)
        
        # list the .readids files in the output directory
        readids_db_fn_l = os.listdir(hitsclip_discardUnclipped_fq_dir())
        readids_db_fn_l = [s for s in readids_db_fn_l if s.split('.')[-1] == 'readids']
        readids_db_fn_l.sort()
        print('Checking that readid db files were created during tests...')
        for fn in readids_db_fn_l:
            print(fn)
        
        expected_fn_l = ['SLBP_CLIP_high_MW_high_MNase.HISTONLY.discardUnclipped.fq.readids',
                         'SLBP_CLIP_high_MW_low_MNase.HISTONLY.discardUnclipped.fq.readids', 
                         'SLBP_CLIP_low_MW_high_MNase.HISTONLY.discardUnclipped.fq.readids', 
                         'SLBP_CLIP_low_MW_low_MNase.HISTONLY.discardUnclipped.fq.readids', 
                         'mock_CLIP_high_MW_high_MNase.HISTONLY.discardUnclipped.fq.readids', 
                         'mock_CLIP_low_MW_high_MNase.HISTONLY.discardUnclipped.fq.readids']
        
        # assert that the expected files have been produced
        self.assertEqual(readids_db_fn_l, expected_fn_l, 
                         'unexpected file output:{}'.format('\n' + \
                         '\n'.join(readids_db_fn_l)))
        
        return

#NOTE: run clippyl.sample_data.paths.remove_test_files to remove these files

if __name__ == '__main__':
    unittest.main()
