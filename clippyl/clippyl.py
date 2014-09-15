import sys
import argparse

from .hitsclip_graphics import hitsclip_graphics_cli
from .build_cleaved_read_db import build_read_db_cli

class Usage(Exception):
    def __init__(self, exitStat):
        self.exitStat = exitStat

#TODO: build the CLI for clippyl ... call the hitclip pipeline routine from here (and other routines).
#http://gehrcke.de/2014/02/distributing-a-python-command-line-application/
#https://github.com/jgehrcke/python-cmdline-bootstrap/blob/master/bootstrap/__main__.py
#https://python-packaging-user-guide.readthedocs.org/en/latest/current.html

#subcommand docs in argparse https://docs.python.org/dev/library/argparse.html#sub-commands

#placeholder

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            d = '''This CLI exposes the full functionality of clippyl'''
            
            #TODO: provide list of possible functions
            #eg hitsclip_graphics, parclip_graphics, iclip_graphics, build_ce_db
            # create the top-level parser
            parser = argparse.ArgumentParser(description = d)
            # subcommand docs in argparse https://docs.python.org/dev/library/argparse.html#sub-commands
            subparsers = parser.add_subparsers()
            #TODO: use generic in and out commands and  inheritance
            
            ####SUBPARSER HITS-CLIP GRAPHICS
            # create the subparser for the "hitsclip_graphics" command
            parser_hcg = subparsers.add_parser('hitsclip_graphics')
            
            #bam_fp_l, required
            parser_hcg.add_argument('bam_files', nargs='+')
            
            #readid_db_fp_l, optional, there must be one cleav_file per bam_file
            parser_hcg.add_argument('--cleav_files', nargs='+')
            
            #query_bed_fp, required
            parser_hcg.add_argument('-q', '--query', nargs='?')
            
            #ciselement_bed_fp_l
            parser_hcg.add_argument('-e', '--ciselements', nargs='+')
            
            #out_pdf_fp #TODO: use pwd as default
            parser_hcg.add_argument('-o', '--output', nargs='?')
            
            parser_hcg.set_defaults(func=hitsclip_graphics_cli)
            
            ####SUBPARSER BUILD CLEAVAGE DB
            # create the subparser for the "build_cleav_db" command
            parser_hcg = subparsers.add_parser('build_cleav_db')
            
            #fastq files of adapter clipped-only reads; required
            parser_hcg.add_argument('fq_files', nargs='+')
            #TODO: output directory argument
            parser_hcg.set_defaults(func=build_read_db_cli)
            
            # parse the args and call whatever function was selected
            args = parser.parse_args()
            print(args) #debugging
            args.func(args)
            
            ####SUBPARSER HITS-CLIP BED-DUMP
            # create the subparser for the "hitsclip_graphics" command
            parser_hcg = subparsers.add_parser('bed-dump')
            
            #bam_fp_l, required
            parser_hcg.add_argument('bam_files', nargs='+')
            
            #readid_db_fp_l, optional, there must be one cleav_file per bam_file
            parser_hcg.add_argument('--cleav_files', nargs='+')
            
            parser_hcg.set_defaults(func=hitsclip_bed-dump_cli)
            
        except SystemExit as exitStat:
            raise Usage(exitStat)
    
    except Usage as err:
        return err.exitStat

