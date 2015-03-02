"""test bedgraph-dump on sample data"""
import os
import unittest

#TODO:check if build_readid will be called first during test discovery
#build_ReadidSQLite_dbs(fq_fp_l)
#TODO: call main instead (see test_coverage_vector for template)
from clippyl.bedgraph_dump import main
from clippyl.sample_data.paths import (hitsclip_fq_fp_l,
                                       hitsclip_cleav_db_l,
                                       hitsclip_bam_fp_l_discardUnclipped,
                                       hitsclip_bam_fp_l,
                                       hitsclip_bam_dir,
                                       hitsclip_n_mapped_reads_l,
                                      )

class BedgraphDumpTestCase(unittest.TestCase):
    
    def runTest(self):
        
        argv = hitsclip_bam_fp_l() + \
               ['--discardUnclipped_db'] + hitsclip_cleav_db_l()
        print(argv) #debugging
        main(argv = argv)
        
        # list the .readids files in the output directory
        bg_fn_l = os.listdir(hitsclip_bam_dir())
        bg_fn_l = [s for s in bg_fn_l if s.split('.')[-1] == 'bg']
        bg_fn_l.sort()
        print('Checking that bedgraph files were created during tests...')
        for fn in bg_fn_l:
            print(fn)
        
        expected_fn_l = ['SLBP_CLIP_high_MW_high_MNase.HISTONLY.bot_std_1D.bg',
                         'SLBP_CLIP_high_MW_high_MNase.HISTONLY.bot_std_clv.bg',
                         'SLBP_CLIP_high_MW_high_MNase.HISTONLY.bot_std_cov.bg',
                         'SLBP_CLIP_high_MW_high_MNase.HISTONLY.top_std_1D.bg',
                         'SLBP_CLIP_high_MW_high_MNase.HISTONLY.top_std_clv.bg',
                         'SLBP_CLIP_high_MW_high_MNase.HISTONLY.top_std_cov.bg',
                         'SLBP_CLIP_high_MW_low_MNase.HISTONLY.bot_std_1D.bg',
                         'SLBP_CLIP_high_MW_low_MNase.HISTONLY.bot_std_clv.bg',
                         'SLBP_CLIP_high_MW_low_MNase.HISTONLY.bot_std_cov.bg',
                         'SLBP_CLIP_high_MW_low_MNase.HISTONLY.top_std_1D.bg',
                         'SLBP_CLIP_high_MW_low_MNase.HISTONLY.top_std_clv.bg',
                         'SLBP_CLIP_high_MW_low_MNase.HISTONLY.top_std_cov.bg',
                         'SLBP_CLIP_low_MW_high_MNase.HISTONLY.bot_std_1D.bg',
                         'SLBP_CLIP_low_MW_high_MNase.HISTONLY.bot_std_clv.bg',
                         'SLBP_CLIP_low_MW_high_MNase.HISTONLY.bot_std_cov.bg',
                         'SLBP_CLIP_low_MW_high_MNase.HISTONLY.top_std_1D.bg',
                         'SLBP_CLIP_low_MW_high_MNase.HISTONLY.top_std_clv.bg',
                         'SLBP_CLIP_low_MW_high_MNase.HISTONLY.top_std_cov.bg',
                         'SLBP_CLIP_low_MW_low_MNase.HISTONLY.bot_std_1D.bg',
                         'SLBP_CLIP_low_MW_low_MNase.HISTONLY.bot_std_clv.bg',
                         'SLBP_CLIP_low_MW_low_MNase.HISTONLY.bot_std_cov.bg',
                         'SLBP_CLIP_low_MW_low_MNase.HISTONLY.top_std_1D.bg',
                         'SLBP_CLIP_low_MW_low_MNase.HISTONLY.top_std_clv.bg',
                         'SLBP_CLIP_low_MW_low_MNase.HISTONLY.top_std_cov.bg',
                         'mock_CLIP_high_MW_high_MNase.HISTONLY.bot_std_1D.bg',
                         'mock_CLIP_high_MW_high_MNase.HISTONLY.bot_std_clv.bg',
                         'mock_CLIP_high_MW_high_MNase.HISTONLY.bot_std_cov.bg',
                         'mock_CLIP_high_MW_high_MNase.HISTONLY.top_std_1D.bg',
                         'mock_CLIP_high_MW_high_MNase.HISTONLY.top_std_clv.bg',
                         'mock_CLIP_high_MW_high_MNase.HISTONLY.top_std_cov.bg',
                         'mock_CLIP_low_MW_high_MNase.HISTONLY.bot_std_1D.bg',
                         'mock_CLIP_low_MW_high_MNase.HISTONLY.bot_std_clv.bg',
                         'mock_CLIP_low_MW_high_MNase.HISTONLY.bot_std_cov.bg',
                         'mock_CLIP_low_MW_high_MNase.HISTONLY.top_std_1D.bg',
                         'mock_CLIP_low_MW_high_MNase.HISTONLY.top_std_clv.bg',
                         'mock_CLIP_low_MW_high_MNase.HISTONLY.top_std_cov.bg',
                        ]
        
        # assert that the expected files have been produced
        self.assertEqual(bg_fn_l, expected_fn_l, 
                         'unexpected file output:{}'.format('\n' + \
                         '\n'.join(bg_fn_l)))
        
        return

#NOTE: run clippyl.sample_data.paths.remove_test_files to remove these files
#from clippyl.sample_data.paths import remove_test_files
#remove_test_files()

if __name__ == '__main__':
    unittest.main()

