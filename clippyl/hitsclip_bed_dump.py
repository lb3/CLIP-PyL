import pysam
import os

from .sqlite_io import ReadidSQLite
from .vector_factory import hitsclip_vectors_2_bg

def hitsclip_bed_dump_cli(args):
    print(args.bam_files,
          args.cleav_files)
    
    hitsclip_bed_dump(args.bam_files,
                      args.cleav_files)
    
    return

def hitsclip_bed_dump(  bam_fp_l,
                        readid_db_fp_l,
                      ):
    
    '''Dump hitsclip stats to bed files'''
    
    for bam_fp, readid_db_fp in zip(bam_fp_l, readid_db_fp_l):
        
        # pass file handles (fh) to the plotting routine
        # 1. load the indexed bam files via pysam filehandles
        #    note: corresponding .bai files are assumed to be in the same directory
        bam_fh = pysam.Samfile(bam_fp, "rb")
        # get the total number of mapped reads per million for each bam file.
        # these values can be used as a normalizing factor
        norm_factor = bam_fh.mapped/1000000
        #TODO: incorporate norm_mode for normalized raw cov
        
        sample_name = os.path.splitext(os.path.basename(bam_fp))
        
        # 2. connect to the readid databases that hold the readids of the adapter 
        #    clipped reads
        #TODO: check for readid file extension and build new database if not found
        cleaved_readid_db_conn = ReadidSQLite(readid_db_fp)
        
        # instantiate bedgraph file handles
        with open(os.path.splitext(bam_fp)[0]+'top_std_cov.bg','w') as bg_fh_top_strand_cov, \
             open(os.path.splitext(bam_fp)[0]+'bot_std_cov.bg','w') as bg_fh_bot_strand_cov, \
             open(os.path.splitext(bam_fp)[0]+'top_std_clv.bg','w') as bg_fh_top_strand_clv, \
             open(os.path.splitext(bam_fp)[0]+'bot_std_clv.bg','w') as bg_fh_bot_strand_clv, \
             open(os.path.splitext(bam_fp)[0]+'top_std_1D.bg','w') as bg_fh_top_strand_1D, \
             open(os.path.splitext(bam_fp)[0]+'bot_std_1D.bg','w') as bg_fh_bot_strand_1D:
            
            f = hitsclip_vectors_2_bg()
            f( bam_fh,
               sample_name,
               bg_fh_top_strand_cov,
               bg_fh_bot_strand_cov,
               bg_fh_top_strand_clv,
               bg_fh_bot_strand_clv,
               bg_fh_top_strand_1D,
               bg_fh_bot_strand_1D,
               cleaved_readid_db_conn = cleaved_readid_db_conn,
               all_adapter_clipped = False,
               uniq_only = True,
               oneD_rate_mode = True,
               oneD_rate_cutoff = 15,
               cleav_rate_mode = False,
               cleav_rate_cutoff = None,
               max_chunk_size = 1000000 )
    
    return
