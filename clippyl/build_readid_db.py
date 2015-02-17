import os
import sys
import argparse
import time

from clippyl.sqlite_io import ReadidSQLite


class Usage(Exception):
    def __init__(self, exitStat):
        self.exitStat = exitStat

def main(argv=None):
    if argv is None:
        argv = sys.argv
    
    try:
        try:
            
            # create the top-level parser
            d = '''This is the CLI for clippyl's readid ''' +\
                '''database builder.\nUse this tool to produce ''' +\
                '''a clippyl-compatible database containing ''' +\
                '''the IDs of\nreads that contained adapter sequence ''' +\
                '''and were clipped.'''
            parser = argparse.ArgumentParser(description=d)
            
            #fastq files of adapter clipped-only reads; required
            parser.add_argument('fq_files', nargs='+')
            
            #output directory
            parser.add_argument('out_dir', nargs='?', const=None, default=None)
            
            args = parser.parse_args()
            #print(args) #debugging
            
            #TODO: allow argparse from argument file
            #https://docs.python.org/3/library/argparse.html#fromfile-prefix-chars
            
            build_cleaved_read_db(args.fq_files)
        
        except SystemExit as exitStat:
            raise Usage(exitStat)
    
    except Usage as err:
        return err.exitStat

def build_ReadidSQLite_dbs(fp_l, out_dir = None):
    """Build a ReadidSQLite database containing the 
    read IDs from a list of fastq files.
    """
    
    print('#######################################')
    print('extracting readids from fastq file')
    
    for in_fp in fp_l:
        
        # using directory where input files are found as default
        # NOTE: DEFAULT FILE INPUT IS GZIP
        if not out_dir:
            out_dir = os.path.dirname(in_fp)
            file_name, file_ext = os.path.splitext(os.path.basename(in_fp))
            out_db_fp = os.path.join(out_dir, file_name + '.readids')
            print('readids will be written to:')
            print(out_db_fp)
        
        start_time = time.time()
        out_db_fh = ReadidSQLite(out_db_fp)
        n = out_db_fh.input_fastq(in_fp)
        elapsed_time = time.time() - start_time
        print()
        print('The amount of time that elapsed during the process was:')
        print('{0:.2f}'.format(round(elapsed_time,2)) + ' seconds')
        print('The number of reads that were processed is:')
        print(str(n))
    
    print('#######################################')
    
    return

