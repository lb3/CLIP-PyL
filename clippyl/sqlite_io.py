# -*- coding: utf-8 -*-
import os
import sqlite3
import time

from clippyl.flatfile_parsing import (FastqReader,
                                      validate_fastq_file,
                                      get_cluster_coords,
                                      Bed6Reader)

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

class ReadidSQLite(SQLiteBase):
    """
    Prior to genome alignment, the CLIP-seq reads must be pre-processed. The 
    pre-processing tasks include quality trimming and sequencing adapter removal.
    Note that the adapter sequence marks the 3' terminus of the RNA fragment of 
    origin. In contrast, the fragments that do not contain adapter sequence merely 
    tag the 5' end of a larger fragment and the 3' terminus of the RNA fragment 
    cannot be known. This software package includes methods to map the 
    RNA fragment termini, which represent cleavage sites that represent either 
    natural cleavages or nuclease cleavages that occur during the controlled 
    nuclease cleavage step of the CLIP-seq library prepartion protocol. Therefore, 
    the list of reads that were adapter trimmed must be accessible.
    
    This module contains a custom SQLite database class that is used to store 
    read_id keys for fast membership testing. The class takes a fastq file 
    containing the set of adapter-clipped reads. This fastq file is produced 
    upstream during pre-processing. Typically, I run the AdapterClipper software 
    from the FastX toolkit with the argument "discardUnclipped" to produce this 
    file. This custom SQLite database class contains methods to import the 
    read ids from the adapter-clipped read file and parse their cluster 
    coordinates, which serve as unique integer keys. This enables the user to 
    align the total set of all reads and then perform cleavage mapping from 
    that bam file by using the SQLite database to identify the adapter-clipped 
    reads.
    """
    
    def input_fastq(self, in_fp):
        """
        Create a SQLite database from the input fastq file.
        """
        #TODO: does not overwrite existing table. you should catch the error and 
        # warn user to delete the file, then quit.
        #TODO: call get_connection
        
        print('...building the database of adapter-clipped reads')
        print(' from the fastq file at:')
        print(in_fp)
        
        # validate the file format and read id format. this is necessary because
        # the read id must be valid to ensure that the cluster coords will be 
        # properly parsed
        validate_fastq_file(in_fp)
        
        with FastqReader(in_fp) as fq_gen:
            
            stat_dict = {'readid_cnt' : 0}
            
            def row_gen():
                for d in fq_gen:
                    t = get_cluster_coords(d['TOPID'])
                    stat_dict['readid_cnt'] += 1
                    if stat_dict['readid_cnt'] % 1000000 == 0:
                        print('{0} reads written'.format(str(stat_dict['readid_cnt'])))
                    yield t
            
            self.c.execute('''CREATE TABLE read_coords(
                                   lane INT not null, 
                                   tile INT not null, 
                                   x_coord INT not null, 
                                   y_coord INT not null,
                                   PRIMARY KEY(lane, tile, x_coord, y_coord))''')
            
            sql_stmnt = '''INSERT INTO read_coords(lane, 
                                                   tile, 
                                                   x_coord, 
                                                   y_coord) VALUES (?,?,?,?)'''
            
            try:
                self.c.executemany(sql_stmnt, row_gen())
            except sqlite3.IntegrityError:
                #TODO: make the original error quit the program too
                print( fq_gen.d['TOPID'] )
                print( fq_gen.get_cluster_coords() )
                print( stat_dict['readid_cnt'] )
        
        self.conn.commit()
        self.conn.close()
        
        return stat_dict['readid_cnt']
    
    def readid_lookup(self, readid):
        '''
        Check if the readid given exists in the database.
        '''
        
        t = get_cluster_coords(readid)
        
        self.c.execute('''SELECT * FROM read_coords
                            WHERE lane = ?
                            AND tile = ?
                            AND x_coord = ?
                            AND y_coord = ?''', t)
        
        return self.c.fetchall()

class Bed6SQLite(SQLiteBase):
    
    def input_bed6(self, fp):
        
        self.c.execute('''CREATE TABLE bed6_intervals
                             (
                                ref TEXT, 
                                start INT, 
                                end INT, 
                                name TEXT, 
                                score INT, 
                                strand TEXT, 
                                PRIMARY KEY(strand, ref, start, end)
                             )''')
        
        stat_dict = {}
        stat_dict['n_of_intervals'] = 0 #total number of intervals
        
        with Bed6Reader(fp) as bed6_gen:
            
            def row_gen():
                
                for d in bed6_gen:
                    
                    stat_dict['n_of_intervals'] += 1
                    
                    yield (d['ref'], 
                           d['start'], 
                           d['end'], 
                           d['name'], 
                           d['score'], 
                           d['strand'])
            
            self.c.executemany('''INSERT INTO bed6_intervals VALUES (?,?,?,?,?,?)''', row_gen())
        
        print('Read', stat_dict['n_of_intervals'], 'intervals')
        self.conn.commit()
        return
    
    def ome_coord_lookup(self, query_coord):
        '''
        Look for all intervals that overlap with the query_coord.
        Return a list of dictionaries, one for each bed6 entry
        that intersects with the query_siteTup.
        '''
        # q_ prefix denotes "query"
        q_ref, q_start, q_end, q_strand = query_coord
        
        t = (q_strand, q_ref, q_start, q_end, q_start, q_end)
        
        self.c.execute('''SELECT * FROM bed6_intervals
                            WHERE strand = ?
                            AND ref = ?
                            AND (start BETWEEN ? AND ? OR end BETWEEN ? AND ?)''', t)
        
        return self.c.fetchall()

