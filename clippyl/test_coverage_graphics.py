import os

from clippyl.build_readid_db import build_ReadidSQLite_dbs
from clippyl.coverage_graphics import hitsclip_graphics
from clippyl.sample_data.paths import (hitsclip_fq_fp_l,
                                       hitsclip_cleav_db_l,
                                       hitsclip_bam_fp_l_discardUnclipped,
                                       hitsclip_bam_fp_l,
                                       hitsclip_n_mapped_reads_l,
                                       histone_gene_bed_fp,
                                       histone_sl_bed_fp,
                                       histone_sl_bed_sl3_fp
                                      )


if __name__ == '__main__':
    
    hitsclip_fq_fp_l()
    hitsclip_cleav_db_l()
    #TODO:check if build_readid will be called first during test discovery
    #build_ReadidSQLite_dbs(fq_fp_l)
    hitsclip_bam_fp_l_discardUnclipped()
    hitsclip_bam_fp_l()
    hitsclip_n_mapped_reads_l()
    histone_gene_bed_fp()
    histone_sl_bed_fp()
    histone_sl_bed_sl3_fp()
    #hitsclip_graphics(

