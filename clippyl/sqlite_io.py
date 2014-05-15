# -*- coding: utf-8 -*-
"""
Prior to genome alignment, the CLIP-seq reads must be pre-processed. The 
pre-processing tasks include quality trimming and sequencing adapter removal.
Note that the adapter sequence marks the 3' terminus of the RNA fragment of 
origin. In contrast, the fragments that do not contain adapter sequence merely 
tag the 5' end of a larger fragment and the 3' terminus of the RNA fragmnet, 
therefore, cannot be known. This software package includes methods to map the 
RNA fragment termini, which represent cleavage sites that represent either 
natural cleavages or nuclease cleavages that occur during the controlled 
nuclease cleavage step of the CLIP-seq library prepartion protocol. Therefore, 
the list of reads that were adapter trimmed must be accessible.

This module contains a custom SQLite database class that is used to store 
read_id keys for fast membership testing. The class takes a fastq file 
containing the set of adapter-clipped reads. This fastq file is produced 
upstream during pre-processing. Typically, I run the AdapterClipper software 
from the FastX toolkit with the argument "discardUnclipped" to produce this 
file. Note that the adapter-clipped reads are a subset of the total set of 
reads. This custom SQLite database class contains methods to import the 
read ids from the adapter-clipped read file and store unique integer keys 
corresponding to the read ids. This enables the user to align the total set 
of all reads and then perform cleavage mapping from that bam file by using the 
SQLite database to identify the adapter-clipped reads.
"""
#TODO: adhere to docstring conventions docs:
#http://sphinxcontrib-napoleon.readthedocs.org/en/latest/example_google.html
#http://legacy.python.org/dev/peps/pep-0257/
#http://legacy.python.org/dev/peps/pep-0008/

import pysam
#docs:
#pysam development branch: https://github.com/pysam-developers/pysam
#pysam docs: http://pysam.readthedocs.org/en/latest/
#pysam docs hosted here too: http://www.cgat.org/~andreas/documentation/pysam/api.html
import sqlite3

class SQLiteBase():
    '''
    A generic sqlite3 database interface.
    '''
    def __init__(self, fp):
        self.fp = fp
        
        self.conn = None
        #TODO: check for file and if it exists, warn of overwrite.
        #(sqlite throws an error if the table already exists)
        self.conn = sqlite3.connect(fp)
        self.c = self.conn.cursor()
        self.c.execute('SELECT SQLITE_VERSION()')
        self.data = self.c.fetchone()
        print('Loaded SQLite database file at:')
        print(fp)
        print('SQLite version: {0}'.format(str((self.data))))
        #TODO: print database structure (tables with # cols and rows)
    
    def __enter__(self):
        return self
    
    def __exit__(self, *exc_info):
        
        self.conn.commit()
        self.conn.close()

def input_fastq(in_fp, out_db_fp):
    """
    Create a SQLite database from the input fastq file.
    """
    #TODO: make this an input cleaved read
    #id method instead and allow input of a list of ids instead
    
    #TODO: overwrite existing file, warn user
    #TODO: call get_connection
    
    print('...building the database of adapter-clipped reads')
    print(' from the fastq file at:')
    print(in_fp)
    
    from clippyl.flatfile_parsing import (FastqReader,
                                          detect_fq_pair_info,
                                          validate_fastq_format)
    
    # validate the file format and read id format. this is necessary because
    # the read id must be valid to ensure that the integer key will be unique
    validate_fastq_format(in_fp)
    
    # bwa seems to remove the matepair tag from the id, therefore it must
    # be removed here to enable comparison downstream
    print('Checking for mate-pair tag in fastq...')
    matepair_tag_bool = False
    if detect_fq_pair_info(in_fp):
        matepair_tag_bool = True
    else:
        pass
    
    fq_gen = FastqReader(in_fp)
    
    stat_dict = {'readid_cnt' : 0}
    
    out_db_fh = SQLiteBase(out_db_fp)
    
    def row_gen():
        for d in fq_gen:
            t = (readid_int_keygen(d['TOPID'], 
                                   strip_matepair_info = matepair_tag_bool),
                )
            stat_dict['readid_cnt'] += 1
            if stat_dict['readid_cnt'] % 1000000 == 0:
                print('{0} reads written'.format(str(stat_dict['readid_cnt'])))
            yield t
    
    out_db_fh.c.execute('''CREATE TABLE adapter_reads(
                           read_id_key INT not null,
                           PRIMARY KEY(read_id_key))''')
    
    sql_stmnt = '''INSERT INTO adapter_reads(read_id_key) VALUES (?)'''
    out_db_fh.c.executemany(sql_stmnt, row_gen())
    
    out_db_fh.conn.commit()
    out_db_fh.conn.close()
    #TODO: test speed of indexes after the database is full instead of
    # declaring the primary key during db creation
    
    return stat_dict['readid_cnt']

def readid_int_keygen(readid, strip_matepair_info = True):
    '''Generates a unique integer from the readid'''
    # Must be a valid read-id format to ensure that the cluster coordinates
    # can be parsed for use downstream during generation of the unique integer
    # key. Must inspect the readid upstream with validate_readid (see below)
    if strip_matepair_info:
        readid
    
    #TODO: strip matepair if found
    # concatenating all integers will create a unique integer key
    readid_key = ''.join([c for c in readid if c.isdigit()])
    
    return readid_key


#TODO: actually code th main function here and link it up to the
# sample data in the clippyl directory structure so that it can be
# run as a test.
##===================MAIN FUNC AREA
import time

print('#######################################')
print('extracting readids from fastq file')

in_fp = '/storage/Ziggy_BigGuy/LB_Bioinformatics/data_Projects/Brooks_HITS_CLIP_SLBP/cleavage_site_mapping/preprocessed_Fastq_discardUnclipped/s_1xS01_sequence.PP.fastq'
#in_fp = '/home/lbthrice/Desktop/fastq_sample/clippedOnly/SRR189782.PP.fastq'
out_db_fp = 'test_qname2_TILEnum2.dat'
start_time = time.time()
n = input_fastq(in_fp, out_db_fp)
elapsed_time = time.time() - start_time

print('qnames were written and indexed to:')
print(out_db_fp)
print('The amount of time that elapsed during the process was:')
print('{0:.2f}'.format(round(elapsed_time,2)) + ' seconds')
print('The number of reads that were processed is:')
print(str(n))
print('#######################################')



