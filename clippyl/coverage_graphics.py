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
#TODO: check for readid file arguments and assume all adpater-clipped if not
###### found. make cleavage plot off by default if not readid file is there

class Usage(Exception):
    def __init__(self, exitStat):
        self.exitStat = exitStat

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            
            # create the top-level parser
            d = '''This is the CLI for generating coverage graphics using clippyl'''
            parser = argparse.ArgumentParser(description=d)
            
            #bam_fp_l, required
            parser.add_argument('-i', '--input_bam_files', nargs='+', required=True)
            
            #query_bed_fp, required
            parser.add_argument('-q','--query_bed_file', nargs='?', required=True)
            
            #stranded
            parser.add_argument('--not_stranded', action='store_true')
            
            #mutually exclusive options allow user to choose a normalization mode
            group = parser.add_mutually_exclusive_group()
            group.add_argument('--n_mapped_reads', nargs='+', type=int)
            group.add_argument('--auto_norm', action='store_true')
            
            #readid_db_fp_l, (optional kwarg)
            #note: there must be one discardUnclipped_db per bam_file
            parser.add_argument('--discardUnclipped_db', nargs='+')
            
            #ciselement_bed_fp_l
            parser.add_argument('-e', '--ciselements', nargs='+')
            
            #out_pdf_fp #TODO: use pwd as default
            parser.add_argument('-o', '--output', nargs='?')
            
            #only hits-clip is currently supported (optional kwarg)
            #TODO: incorporate code for par-clip and iclip
            parser.add_argument('--clipseq_method', choices=['hits-clip',],
                                default='hits-clip')
            
            # parse the args and call whatever function was selected
            args = parser.parse_args(argv)
            
            if args.clipseq_method == 'hits-clip':
                hitsclip_graphics(args)
            else:
                #TODO: incorporate par-clip and iclip options
                pass
            
            print(args) #debugging
        
        except SystemExit as exitStat:
            raise Usage(exitStat)
    
    except Usage as err:
        return err.exitStat

def hitsclip_graphics( args ):
    '''Batch graphics routine for hitsclip data'''
    
    print(args)
    
    bam_fp_l = args.input_bam_files
    not_stranded = args.not_stranded
    #TODO: insert n_mapped_reads attribute and check that None state is handled
    n_mapped_reads= args.n_mapped_reads
    auto_norm = args.auto_norm
    readid_db_fp_l = args.discardUnclipped_db
    query_bed_fp = args.query_bed_file
    ciselement_bed_fp_l = args.ciselements
    out_pdf_fp = args.output
    
    if not out_pdf_fp:
        #print('pwd', os.getcwd()) #debugging
        out_pdf_fp = os.path.join(os.getcwd(), 'CLIP-PyL_coverage_graphics.pdf')
    
    if not_stranded:
        stranded = False
        print('PROCESSING AS NON-STRANDED AKA STRAND-AGNOSTIC LIBRARY')
    else:
        stranded = True
    
    # parse data labels from input bam file names
    def parse_label(filename):
        return filename.split('.')[0]
    
    print('bamfpl')
    print(bam_fp_l)
    label_l = [parse_label(os.path.basename(fp)) for fp in bam_fp_l]
    
    # pass file handles (fh) to the plotting routine
    # 1. load the indexed bam files via pysam filehandles
    #    note: corresponding .bai files are assumed to be in the same directory
    bam_fh_l = [pysam.Samfile(fp, "rb") for fp in bam_fp_l]
    
    # OPTIONAL
    # mutually exclusive with options for normaliztion
    if auto_norm:
        # Get the total number of mapped reads per million for each bam file
        # and use those values as the normalizing factors (the read counts 
        # will be divided by million mapped reads in each file
        norm_factor_l = [bam_fh.mapped/1000000  for bam_fh in bam_fh_l]
    elif n_mapped_reads:
            for i in n_mapped_reads: #debugging
                assert type(i) == int #debugging
            norm_factor_l = [int(i)/1000000 for i in n_mapped_reads]
    else:
        norm_factor_l = None
    
    # OPTIONAL
    # connect to the readid databases that hold the readids of the adapter 
    # clipped reads. this is an optional kwarg upstream at the main CLI
    if readid_db_fp_l:
        readid_db_fh_l = [ReadidSQLite(fp) for fp in readid_db_fp_l]
    else:
        readid_db_fh_l = None
    
    # load the bed file containing query intervals into a generator function
    bed_gen = Bed6Reader(query_bed_fp)
    
    # OPTIONAL
    # load the cis element interval database (sqlite3 file format)
    # check if file has .sl3 file extension, if not then build sl3 db
    # and connect to it
    if ciselement_bed_fp_l:
        ciselement_db_fh_l = []
        for fp in ciselement_bed_fp_l:
            if fp.split('.')[-1] == 'sl3':
                ciselement_db_fh_l.append(Bed6SQLite(fp))
            else:
                print('building ciselement db at:')
                print(fp)
                fh = Bed6SQLite(fp + '.sl3')
                fh.input_bed6(fp)
                ciselement_db_fh_l.append(fh)
        # collect user-defined cis element labels from params file
        #TODO: parse this from input file names
        ciselement_label_l = ['histone_SL', ]
    else:
        ciselement_db_fh_l = None
        ciselement_label_l = None
    
    # load the graphics output file handle where plots will be aggregated
    # docs for multi-page pdf output found here: 
    # http://matplotlib.sourceforge.net/faq/howto_faq.html#save-multiple-plots-to-one-pdf-file
    pp = PdfPages(out_pdf_fp)
    
    for i in bed_gen:
        
        #debugging
        print( bed_gen,
               bam_fh_l,
               readid_db_fh_l,
               label_l,
               norm_factor_l,
               ciselement_db_fh_l,
               ciselement_label_l )
        
        fig = hits_clip_plot(   bed_gen,
                                bam_fh_l,
                                readid_db_fh_l = readid_db_fh_l,
                                label_l = label_l,
                                norm_factor_l = norm_factor_l,
                                ciselement_db_fh_l = ciselement_db_fh_l,
                                ciselement_label_l = ciselement_label_l,
                                ciselement_color_l = None,
                                flank = 100,
                                stranded = stranded )
        
        pp.savefig()
    pp.close()
    
    return

