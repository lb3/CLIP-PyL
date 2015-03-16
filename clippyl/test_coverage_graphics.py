'''test coverage graphics program'''
import os
import unittest

from clippyl.build_readid_db import build_ReadidSQLite_dbs
from clippyl.coverage_graphics import main
from clippyl.sample_data.paths import (hitsclip_fq_fp_l,
                                       hitsclip_cleav_db_l,
                                       hitsclip_bam_fp_l_discardUnclipped,
                                       hitsclip_bam_fp_l,
                                       hitsclip_n_mapped_reads_l,
                                       histone_gene_bed_fp,
                                       histone_sl_bed_fp,
                                       histone_sl_bed_sl3_fp
                                      )

class Coverage_Graphics_Basic_TC(unittest.TestCase):
    
    def test_coverage_graphics(self):
        # generate coverage graphics for adapter-clipped reads
        argv=['-i',] + hitsclip_bam_fp_l_discardUnclipped() + \
             ['--n_mapped_reads'] + hitsclip_n_mapped_reads_l(out_type='string') + \
             ['-q', histone_gene_bed_fp(), 
              '--output', 'CLIP-PyL_graphics_test.pdf']
        #print('argv:', argv) #debugging
        main(argv=argv)
        
        return

class CoverageGraphicsTestCase(unittest.TestCase):
    
    def test_coverage_graphics(self):
        # generate coverage graphics for adapter-clipped reads
        argv=['-i',] + hitsclip_bam_fp_l_discardUnclipped() + \
             ['--n_mapped_reads'] + hitsclip_n_mapped_reads_l(out_type='string') + \
             ['-q', histone_gene_bed_fp(), 
              '-e', histone_sl_bed_fp(), 
              '--output', 'CLIP-PyL_graphics_test_1.pdf']
        #print('argv:', argv) #debugging
        main(argv=argv)
        
        return
    
    def test_graphics_w_readid_db(self):
        # generate coverage graphics that include all reads
        argv=['-i',] + hitsclip_bam_fp_l() + \
             ['--auto_norm'] + \
             ['--discardUnclipped_db',] + hitsclip_cleav_db_l() + \
             ['-q', histone_gene_bed_fp(), 
              '-e', histone_sl_bed_sl3_fp(), 
              '--output', 'CLIP-PyL_graphics_test_2.pdf']
        #print('argv:', argv) #debugging
        main(argv=argv)
        
        return

#NOTE: You can run clippyl.sample_data.paths.remove_test_files 
#      to remove these test files; template code:
#from clippyl.sample_data.paths import remove_test_files
#remove_test_files()

if __name__ == '__main__':
    unittest.main()

