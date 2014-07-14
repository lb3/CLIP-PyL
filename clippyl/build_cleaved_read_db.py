import time
import os

from .sqlite_io import ReadidSQLite
from .flatfile_parsing import FastqReader

def build_read_db_cli(args):
    '''this is a function to map the namespace object from the argparse CLI'''
    build_cleaved_read_db(args.fq_files)
    return

def build_cleaved_read_db(fp_l):
    '''build a sqlite database containing the coordinates of the cleaved reads 
       in the fastq files'''
    
    print('#######################################')
    print('extracting readids from fastq file')
    
    for in_fp in fp_l:
        
        #determine output file name for the sqlite3 db
        # parameterize out_dir. use 
        #TODO: use directory where inputy files are found as default instead
        out_dir = os.getcwd()
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

#print('#######################################')
#print('testing lookup function')
#from clippyl.sqlite_io import ReadidSQLite
#from clippyl.flatfile_parsing import FastqReader

#fq_fp = '/home/lbthrice/Desktop/fastq_sample/clippedOnly/s_1xS01_sequence.PP.fastq'
#db_fp = './s_1xS01_sequence.PP.readids'

#with FastqReader(fq_fp) as fq_gen, ReadidSQLite(db_fp) as readid_db, open('t.txt', 'w') as out_fh:
#    for d in fq_gen:
#        out_fh.write( repr(readid_db.readid_lookup(d['TOPID'])) + '\n')
#print('#######################################')
