import pysam
import os

from .flatfile_parsing import Bed6Reader
from .sqlite_io import ReadidSQLite

def hitsclip_bed-dump_cli(args):
    print(args.bam_files,
          args.cleav_files,
          args.output)
    
    hitsclip_bed-dump(args.bam_files,
                      args.cleav_files,
                      args.output)
    return

def hitsclip_bed-dump(  bam_fp_l,
                        readid_db_fp_l,
                        out_pdf_fp ):
    
    '''Dump hitsclip stats to bed files'''
    
    # parse data labels from input bam file names
    def parse_label(filename):
        return filename.split('.')[0]
    
    label_l = [parse_label(os.path.basename(fp)) for fp in bam_fp_l]
    
    # pass file handles (fh) to the plotting routine
    # 1. load the indexed bam files via pysam filehandles
    #    note: corresponding .bai files are assumed to be in the same directory
    bam_fh_l = [pysam.Samfile(fp, "rb") for fp in bam_fp_l]
    
    # get the total number of mapped reads per million for each bam file.
    # these values can be used as a normalizing factor
    norm_factor_l = [bam_fh.mapped/1000000  for bam_fh in bam_fh_l]
    
    # 2. connect to the readid databases that hold the readids of the adapter 
    #    clipped reads
    #TODO: check for readid file extension and build new database if not found
    readid_db_fh_l = [ReadidSQLite(fp) for fp in readid_db_fp_l]
    
    # 3. load the graphics output file handle where plots will be aggregated
    # docs for multi-page pdf output: http://matplotlib.sourceforge.net/faq/howto_faq.html#save-multiple-plots-to-one-pdf-file
    
    for i in bed_gen:
        
        fig = hits_clip_plot(   bed_gen,
                                bam_fh_l,
                                readid_db_fh_l = readid_db_fh_l,
                                norm_factor_l = norm_factor_l,
                                stranded = True )
        
        pp.savefig()
    pp.close()
    return
