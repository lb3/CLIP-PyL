import pysam
import os

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from .flatfile_parsing import Bed6Reader
from .sqlite_io import ReadidSQLite, Bed6SQLite
from .vector_factory import build_hitsclip_vectors
from .mpl_graphics import hits_clip_plot

def hitsclip_graphics_cli(args):
    print(args.bam_files,
          args.cleav_files,
          args.query,
          args.ciselements,
          args.output)
    
    hitsclip_graphics(args.bam_files,
                      args.cleav_files,
                      args.query,
                      args.ciselements,
                      args.output)
    return

def hitsclip_graphics(  bam_fp_l,
                        readid_db_fp_l,
                        query_bed_fp,
                        ciselement_bed_fp_l,
                        out_pdf_fp ):
    '''Batch graphics routine for hitsclip data'''
    
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
    
    # 3. load the bed file containing query intervals into a generator function
    bed_gen = Bed6Reader(query_bed_fp)
    
    # 4. load the cis element interval database (sqlite3 file format)
    #    check if file has .sl3 file extension, if not then build sl3 db
    #    and connect to it
    ciselement_db_fh_l = []
    for fp in ciselement_bed_fp_l:
        if fp.split('.')[-1] == 'sl3':
            ciselement_db_fh_l.append(Bed6SQLite(fp))
        else:
            print('building ciselement db')
            print(fp)
            fh = Bed6SQLite(fp + '.sl3')
            fh.input_bed6(fp)
            ciselement_db_fh_l.append(fh)
    # collect user-defined cis element labels from params file
    #TODO: parse this from input file names
    ciselement_label_l = ['histone_SL', ]
    
    # 5. load the graphics output file handle where plots will be aggregated
    # docs for multi-page pdf output: http://matplotlib.sourceforge.net/faq/howto_faq.html#save-multiple-plots-to-one-pdf-file
    pp = PdfPages(out_pdf_fp)
    
    for i in bed_gen:
        
        fig = hits_clip_plot(   bed_gen,
                                bam_fh_l,
                                readid_db_fh_l = readid_db_fh_l,
                                label_l = label_l,
                                norm_factor_l = norm_factor_l,
                                ciselement_db_fh_l = ciselement_db_fh_l,
                                ciselement_label_l = ciselement_label_l,
                                ciselement_color_l = None,
                                flank = 100,
                                stranded = True )
        
        pp.savefig()
    pp.close()
    return
