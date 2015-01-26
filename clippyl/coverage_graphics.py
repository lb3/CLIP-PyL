import os
import sys
import argparse

import pysam

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from clippyl.flatfile_parsing import Bed6Reader
from clippyl.sqlite_io import ReadidSQLite, Bed6SQLite
from clippyl.mpl_graphics import hits_clip_plot

#TODO: emulate bedgraph-dump behavior and make it default assume all adaper clipped
#eg fix this (copied from below)
    # 2. connect to the readid databases that hold the readids of the adapter 
    #    clipped reads
    #TODO: check for readid file extension and build new database if not found

class Usage(Exception):
    def __init__(self, exitStat):
        self.exitStat = exitStat

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            
            #TODO: provide list of possible functions
            # create the top-level parser
            d = '''This is the CLI for generating coverage graphics using clippyl'''
            parser = argparse.ArgumentParser(description=d)
            
            #bam_fp_l, required
            parser.add_argument('bam_files', nargs='+')
            
            #query_bed_fp, required
            parser.add_argument('-q', '--query', nargs='?')
            
            #readid_db_fp_l, (optional kwarg)
            #note: there must be one cleav_file per bam_file
            parser.add_argument('--cleav_files', nargs='+')
            
            #ciselement_bed_fp_l
            parser.add_argument('-e', '--ciselements', nargs='+')
            
            #out_pdf_fp #TODO: use pwd as default
            parser.add_argument('-o', '--output', nargs='?')
            
            #only hits-clip is currently supported (optional kwarg)
            #TODO: incorporate code for par-clip and iclip
            parser.add_argument('--clipseq_method', choices=['hits-clip',],
                                default='hits-clip')
            
            # parse the args and call whatever function was selected
            args = parser.parse_args()
            
            if args.clipseq_method == 'hits-clip':
                hitsclip_graphics(args.bam_files,
                                  args.cleav_files,
                                  args.query,
                                  args.ciselements,
                                  args.output)
            else:
                #TODO: incorporate par-clip and iclip options
                pass
            
        except SystemExit as exitStat:
            raise Usage(exitStat)
    
    except Usage as err:
        return err.exitStat

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

