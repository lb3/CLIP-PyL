#!/usr/bin/python3
# -*- coding: utf-8 -*-

#cims_db.py
#cleavage_db.py
#raw_coverage_db.py

#NOTE: this function ONLY includes uniquely aligning reads
# this function computes raw coverage
#TODO: add SamGen.isMapped() logic and totalMappedReads dict key as 
#shown in bowtieSam2rawCoverCD
#TODO: detect sam file type with SamReader method, invoke proper parse scheme

import struct

class PositionalShelf():
    """
    API DEVELOPMENT WARNING: ALWAYS EXPLICITLY DECLARE BOTH STRAND KEYS BECAUSE
    SOME FUNCTIONS INJECT A library_stats KEY AT THE TOP LEVEL; DO NOT SIMPLY 
    ITERATE OVER THE TOP LEVEL KEYS OR YOU WILL GET A TYPE-ERROR.
    a dictionary-based storage structure for storing data about chromosome-positions
    the input object is always a tuple with the following data type:
    (str(reference), int(startCoord), int(endCoord), str(strand), dataobject)
    the reference string is most often the chromosome but could be gene model name name or something else
    the coordinates expected are UCSC style (right-open zero based; (startCoord, endCoord])
    add data to the dictionary by providing an input object to the addSiteTup method
    the dataobject data type is defined in the dataType argument of the constructor
    the data stored in the dataObject slot is integer by default therefore arithematic operators occur 
    """
    
    def __init__(self, fp, wb_bool = True, dataType = int, write_chunk = 10000):
        
        import shelve
        #caution, if this is an existing db then you must pass the 
        #filepath stripped of it's extension so that the shelve 
        #detects the db-type. if it is a new file then the file extension
        #will be appended automatically
        self.d = shelve.open(fp, writeback = wb_bool)
        
        #detect if the dictionary is occupied
        if '+' in self.d:
            #TODO: write a dataType key and read it upon opening to
            #set the datatype param properly
            pass
        else:
            self.d['+'] = {}
            self.d['-'] = {}
        
        #TODO: this is kind of silly, as datatype class can be realized at addSiteTup
        self.dataType = dataType
        
        #instantiate counter for sync timing
        self.write_count = 0
        self.write_chunk = write_chunk
        
        return
    
    def __enter__(self):
       return self
    
    def __exit__(self, *exc_info):
       #print('executing data sync for positional db') #debugging
       self.d.sync()
       self.d.close()
    
    def addSiteTup(self, siteTup):
        """
        a method for adding data from a tuple
        by default the addition operation is used for dataType int or float
        by default a dataType of str causes the string objects to be 
        aggregated into lists at each site
        """
        
        reference, startCoord, endCoord, strand, dataObject = siteTup
        
        #perform the correct data operation for the type given
        if self.dataType == int or self.dataType == float:
            #check that the site has been written to the chromDict in self.d storage object
            if reference in self.d[strand]:
                for coord in range(startCoord, endCoord):
                    if coord in self.d[strand][reference]:
                        self.d[strand][reference][coord] += dataObject
                    else:
                        self.d[strand][reference][coord] = dataObject
            else:
                self.d[strand][reference] = {}
                for coord in range(startCoord, endCoord):
                    self.d[strand][reference][coord] = dataObject
        else:
            print('unsupported dataObject dataType!')
            raise IOError
        
        self.write_count += 1
        if self.write_count % self.write_chunk == 0:
            #print('shelf sync event') #debugging
            self.d.sync()
            #print('{0} sites written'.format(str(self.write_count)))
        
        return
    
    def key_exists( self, ref, strand, coord):
        
        if ref in self.d[strand]:
            if coord in self.d[strand][ref]:
                return True
            else:
                return False
        else:
            return False

class SQLiteBase():
    
    def __init__(self, fp):
        
        import sqlite3
        
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


#TODO: build a class method to build a table to hold the 1D offset index 
# and data for downstream permutation
class HitsClipSQLite(SQLiteBase):
    """
    SQLite adapters for HitsClip data.
    """
    
    def input_sam(self, fp, uniq_only = False):
        """
        Input a sam file.
        """
        #NOTE: If you filter for uniq here you have a much smaller db and you save time.
        #HOWEVER! if you are going to do a permutation routine downstream then 
        #      you want to set uniq_only to false and use the uniq_only selection
        #      downstream (in build_basewise_signature) so that all the 1D events 
        #      are available in the sam sqlite for permutations
        
        #TODO: include the function to parse MAXOCC flag in bwa alignments
        #      that were invoked with the maxOcc option defined
        from clipPyL.flatfile_parsing import SamReader
        
        self.c.execute('''CREATE TABLE sam
                     (
                        qname TEXT, 
                        ref TEXT, 
                        start_coord INT, 
                        end_coord INT, 
                        strand TEXT, 
                        cigar TEXT, 
                        uniq_flag INT,
                        PRIMARY KEY(qname)
                     )''')
        
        stat_dict = {}
        stat_dict['n_of_reads'] = 0 #total number of reads is gleaned from sam file
        stat_dict['n_of_mapped_reads'] = 0 #gleaned from sam flagbits
        stat_dict['n_of_unique_aligns'] = 0 #gleaned from sam XT opt field
        
        with SamReader(fp) as sam_gen:
            
            def row_gen():
                
                for d in sam_gen:
                    stat_dict['n_of_reads'] += 1
                    
                    if sam_gen.isMapped():
                        stat_dict['n_of_mapped_reads'] += 1
                        
                        #identify uniquely aligned reads by accessing 
                        #the sam optdict
                        if sam_gen.unique_bool():
                            stat_dict['n_of_unique_aligns'] += 1
                            uniq_flag = 1
                        elif uniq_only:
                            continue
                        else:
                            uniq_flag = None
                        
                        #raw coverage interval
                        reference, startCoord, endCoord, strand = sam_gen.calc_ucsc_interval()
                        #TODO: store 1D offsets for permutation downstream,
                        
                        yield (d['QNAME'], reference, startCoord, endCoord, strand, d['CIGAR'], uniq_flag)
            
            self.c.executemany('''INSERT INTO sam VALUES (?,?,?,?,?,?,?)''', row_gen())
        
        self.conn.commit()
        
        return stat_dict
    
    def input_fastq(self, fp, col = ('id') ):
        """
        Create a SQLite database from the input fastq file.
        Currently, only the id column is supported.
        """
        #TODO: make this an input cleaved read
        #id method instead and allow input of a list of ids instead
        
        #TODO: overwrite existing file, warn user
        #TODO: if useful, implement col = ('id','qual','seq')) 
        #      by putting qual and seq in shelf values.
        #TODO: call get_connection
        
        print('...building the database of adapter-clipped reads')
        print(' from the fastq file at:')
        print(fp)
        print('...building {0} column(s)...'.format( col ))
        
        from clipPyL.flatfile_parsing import FastqReader, detect_fq_pair_info
        
        # bwa seems to remove the matepair tag from the id, therefore it must
        # be removed here to enable comparison
        print('Checking for mate-pair tag fastq...')
        rm_matepair_tag = False
        if detect_fq_pair_info(fp):
            rm_matepair_tag = True
        else:
            pass
        
        fq_gen = FastqReader(fp)
        
        stat_dict = {'cleav_read_cnt' : 0}
        
        
        if col == ('id'):
            
            if rm_matepair_tag:
                def row_gen():
                    for d in fq_gen:
                        t = (d['TOPID'][1:-2], ) #strips leading '@' AND the matpair tag
                        stat_dict['cleav_read_cnt'] += 1
                        if stat_dict['cleav_read_cnt'] % 1000000 == 0:
                            print('{0} reads written'.format(str(stat_dict['cleav_read_cnt'])))
                        yield t
            else:
                def row_gen():
                    for d in fq_gen:
                        t = (d['TOPID'][1:], ) #strips leading '@'
                        stat_dict['cleav_read_cnt'] += 1
                        if stat_dict['cleav_read_cnt'] % 1000000 == 0:
                            print('{0} reads written'.format(str(stat_dict['cleav_read_cnt'])))
                        yield t
            
            self.c.execute('''CREATE TABLE adpt_clipped_reads (
                              read_id TEXT not null, 
                              PRIMARY KEY(read_id))''')
            
            self.c.executemany('''INSERT INTO adpt_clipped_reads VALUES 
                                   (?)''', row_gen())
            
            self.conn.commit()
            #TODO: test speed of indexes after the database is full instead of
            # declaring the primary key during db creation
        
        else:
            raise IOError('clipPyL does not yet support that column structure')
        
        return stat_dict['cleav_read_cnt']
    
    def build_basewise_signature(self, sam_db_fp, cleav_db_fp = None, uniq_only = False):
        '''
        Generate basewise counts for 1D, cleavage and raw coverage from the 
        set of alignements at the bwa-samse output file at sam_db_fp.
        The build process first parses the sam file and exports the calculations
        to positional shelve db files at pos_shelve_fp. The export to shelf 
        is necessary because building the positional sqlite database takes too 
        much time due to the need for UPDATE calls to read each new site. 
        Therefore, the PositonalShelves are written and then copied 
        to the connected sqlite database.
        
        Note: if the cleav_db_fp is set to None then it is assumed that all
        the reads in the sam file are adapater-clipped
        '''
        #TODO: positional db can be a temporary file
        #TODO: sqlite file with sam and fq tables can be temporary
        # the way to do this will be to use temp directory file paths
        # in context managers to generate fp arguments that are passed to 
        # these functions
        #TODO: add an option to input a simple list of clipped reads upstream 
        # at "input_fastq" method. ...which should be changed 
        # to "build_clipped_readid_db" method
        import os
        import tempfile
        import re
        
        sam_db_fn_prefix = os.path.basename(sam_db_fp).split('.')[0]
        temp_dir = tempfile.TemporaryDirectory( suffix=sam_db_fn_prefix, 
                                                dir=os.path.dirname(self.fp) )
        
        CIGARregexP = re.compile(r'([0-9]|[0-9][0-9])([MIDNSHP=X])')
        
        def get1DsiteList(ref, start_coord, strand, cigar_string):
            '''
            A function to glean the list of site tuples comprised of all
            1D events described by the zero-based start_coord argument and
            the cigar string that was copied from the sam file and is
            therefore one-based.
            '''
            oneDSitesList = []
            m = CIGARregexP.findall(cigar_string)
            for idx, tup in enumerate(m):
                posCount, operation = tup
                if operation == 'D':
                    offset = sum([int(t[0]) for t in m[:idx + 1] if t[1] != 'I'])
                    offset -= 1 # convert to zero-based coord
                    oneDsiteTup = (ref, 
                                   start_coord + offset, 
                                   start_coord + offset + 1, 
                                   strand)
                    oneDSitesList.append(oneDsiteTup)
            return oneDSitesList
        
        stat_dict = {}
        stat_dict['n_of_nts_covered'] = 0 # total of nt sites with coverage
        stat_dict['n_of_oneD_sites'] = 0 # total number of nts with oneD
        stat_dict['n_of_cleavage_sites'] = 0 #total number of nts with cleavage
        stat_dict['n_of_mapped_cleav_reads'] = 0 #note: impacted by uniq_only
        stat_dict['n_of_mapped_oneD_reads'] = 0 #note: impacted by uniq_only
        
        if cleav_db_fp:
            
            print('Selecting adapter-clipped reads')
            self.c.execute('''ATTACH ? AS adpt_clipped_reads''',(cleav_db_fp,))
            self.c.execute('''ATTACH ? AS sam''',(sam_db_fp,))
            #select cleaved reads
            if uniq_only:
                # to select reads with non-null uniq_flag here
                self.c.execute('''SELECT ref TEXT, 
                                         start_coord INT, 
                                         end_coord INT, 
                                         strand TEXT, 
                                         cigar TEXT 
                                  FROM sam, adpt_clipped_reads 
                                  WHERE sam.qname = adpt_clipped_reads.read_id 
                                  AND sam.uniq_flag NOT NULL''')
            else:
                self.c.execute('''SELECT ref TEXT, 
                                         start_coord INT, 
                                         end_coord INT, 
                                         strand TEXT, 
                                         cigar TEXT 
                                  FROM sam, adpt_clipped_reads 
                                  WHERE sam.qname = adpt_clipped_reads.read_id''')
            
            print('Tallying adapter-clipped reads')
            # subdirectory and explicity delete it in the script
            with PositionalShelf(os.path.join(temp_dir.name, sam_db_fn_prefix) + '.raw') as p_db_raw, \
                 PositionalShelf(os.path.join(temp_dir.name, sam_db_fn_prefix) + '.cleav') as p_db_cleav, \
                 PositionalShelf(os.path.join(temp_dir.name, sam_db_fn_prefix) + '.oneD') as p_db_oneD:
                
                l = self.c.fetchmany()
                while l:
                    for row in l:
                        stat_dict['n_of_mapped_cleav_reads'] += 1
                        if stat_dict['n_of_mapped_cleav_reads'] % 1000000 == 0:
                            print(stat_dict['n_of_mapped_cleav_reads'], '''alignments of adapter-clipped reads have been processed''')
                        
                        ref, start_coord, end_coord, strand, cigar = row
                        p_db_raw.addSiteTup( (ref, start_coord, end_coord, strand, 1) )
                        
                        #add left cleavage site
                        p_db_cleav.addSiteTup( (ref, start_coord, start_coord + 1, strand, 1) )
                        #add right cleavage site
                        p_db_cleav.addSiteTup( (ref, end_coord - 1, end_coord, strand, 1) )
                        
                        if '1D' in cigar:
                            oneDsite_l = get1DsiteList(ref, start_coord, strand, cigar)
                            if oneDsite_l:
                                stat_dict['n_of_mapped_oneD_reads'] += 1
                                for oneDsite in oneDsite_l:
                                    ref, start_coord, end_coord, strand = oneDsite
                                    p_db_oneD.addSiteTup( (ref, start_coord, end_coord, strand, 1) )
                    
                    l = self.c.fetchmany()
                
                print(stat_dict['n_of_mapped_cleav_reads'], '''alignments of adapter-clipped reads have been processed''')
                
                print('Selecting unclipped reads')
                if uniq_only:
                    self.c.execute('''SELECT ref TEXT, 
                                             start_coord INT, 
                                             end_coord INT, 
                                             strand TEXT, 
                                             cigar TEXT 
                                      FROM sam
                                      LEFT JOIN adpt_clipped_reads ON adpt_clipped_reads.read_id = sam.qname
                                      WHERE adpt_clipped_reads.read_id IS NULL
                                      AND sam.uniq_flag NOT NULL''')
                else:
                    self.c.execute('''SELECT ref TEXT, 
                                             start_coord INT, 
                                             end_coord INT, 
                                             strand TEXT, 
                                             cigar TEXT 
                                      FROM sam
                                      LEFT JOIN adpt_clipped_reads ON adpt_clipped_reads.read_id = sam.qname
                                      WHERE adpt_clipped_reads.read_id IS NULL''')
                
                print('Tallying unclipped reads.')
                n_of_alignments = 0
                l = self.c.fetchmany()
                while l:
                    for row in l:
                        n_of_alignments += 1
                        if n_of_alignments % 1000000 == 0:
                            print(n_of_alignments, '''alignments of unclipped reads have been processed''')
                        
                        ref, start_coord, end_coord, strand, cigar = row
                        p_db_raw.addSiteTup( (ref, start_coord, end_coord, strand, 1) )
                        
                        if '1D' in cigar:
                            oneDsite_l = get1DsiteList(ref, start_coord, strand, cigar)
                            if oneDsite_l:
                                stat_dict['n_of_mapped_oneD_reads'] += 1
                                for oneDsite in oneDsite_l:
                                    ref, start_coord, end_coord, strand = oneDsite
                                    p_db_oneD.addSiteTup( (ref, start_coord, end_coord, strand, 1) )
                    l = self.c.fetchmany()
                
                print(n_of_alignments, '''alignments of unclipped reads have been processed''')
                
                self.c.execute('''DETACH adpt_clipped_reads''')
                self.c.execute('''DETACH sam''')
                
                #perform data sync before write-out
                p_db_raw.d.sync()
                p_db_cleav.d.sync()
                p_db_oneD.d.sync()
                
                #TODO: optimize order of PRIMARY KEY indices
                self.c.execute('''CREATE TABLE clip_signature
                                     (
                                        ref TEXT, 
                                        coord INT, 
                                        strand TEXT, 
                                        raw_coverage INT, 
                                        cleavage INT, 
                                        single_nt_del INT, 
                                        PRIMARY KEY(strand, ref, coord)
                                     )''')
                
                # iterate over 1D sites first, lookup sites in cleav and raw
                # positional db and delete once they are written to the hits-clip 
                # sqlite to maintain unique site index across the shelves.
                def row_gen():
                    print('injecting clip signature values for 1D sites')
                    cnt = 0
                    for strand in p_db_oneD.d.keys():
                        print('loading sites from the {0} strand'.format(strand))
                        for ref in p_db_oneD.d[strand].keys():
                            for coord in p_db_oneD.d[strand][ref].keys():
                                
                                oneD_coverage = p_db_oneD.d[strand][ref][coord]
                                stat_dict['n_of_oneD_sites'] += 1
                                
                                # pull raw coverage value
                                raw_coverage = p_db_raw.d[strand][ref][coord]
                                del p_db_raw.d[strand][ref][coord]
                                stat_dict['n_of_nts_covered'] += 1
                                
                                # pull cleav coverage value
                                if p_db_cleav.key_exists( ref, strand, coord):
                                    cleav_coverage = p_db_cleav.d[strand][ref][coord]
                                    del p_db_cleav.d[strand][ref][coord]
                                    stat_dict['n_of_cleavage_sites'] += 1
                                else:
                                    cleav_coverage = 0
                                
                                #construct the row tuple
                                t = (ref, coord, strand, raw_coverage, cleav_coverage, oneD_coverage)
                                cnt += 1
                                
                                # periodically sync data to writeback key removals
                                if cnt % 1000000 == 0:
                                    p_db_raw.d.sync()
                                    p_db_cleav.d.sync()
                                    print('{0} writes to sqlite'.format(str(cnt)))
                                
                                yield t
                
                # inject sites into sqlite
                self.c.executemany('''INSERT INTO clip_signature VALUES
                                       (?,?,?,?,?,?)''', row_gen())
                
                # sync positional shelf data to writeback key removals
                p_db_raw.d.sync()
                p_db_cleav.d.sync()
                
                # iterate over the remaining cleaved sites, lookup sites in raw
                # positional db and delete once they are written to the hits-clip 
                # sqlite to maintain unique site index across the shelves.
                def row_gen():
                    print('injecting clip signature values for remaining cleavage sites')
                    cnt = 0
                    for strand in p_db_cleav.d.keys():
                        print('loading sites from the {0} strand'.format(strand))
                        for ref in p_db_cleav.d[strand].keys():
                            for coord in p_db_cleav.d[strand][ref].keys():
                                
                                cleav_coverage = p_db_cleav.d[strand][ref][coord]
                                stat_dict['n_of_cleavage_sites'] += 1
                                
                                # pull raw coverage value
                                raw_coverage = p_db_raw.d[strand][ref][coord]
                                del p_db_raw.d[strand][ref][coord]
                                stat_dict['n_of_nts_covered'] += 1
                                
                                #construct the row tuple
                                t = (ref, coord, strand, raw_coverage, cleav_coverage, 0)
                                cnt += 1
                                
                                # periodically sync data to writeback key removals
                                if cnt % 1000000 == 0:
                                    p_db_raw.d.sync()
                                    print('{0} writes to sqlite'.format(str(cnt)))
                                
                                yield t
                
                # inject sites into sqlite
                self.c.executemany('''INSERT INTO clip_signature VALUES
                                       (?,?,?,?,?,?)''', row_gen())
                
                # sync positional shelf data to writeback key removals
                p_db_raw.d.sync()
                
                # iterate over the remaining sites
                def row_gen():
                    print('injecting raw coverage values for remaining sites')
                    cnt = 0
                    for strand in p_db_raw.d.keys():
                        print('loading sites from the {0} strand'.format(strand))
                        for ref in p_db_raw.d[strand].keys():
                            for coord in p_db_raw.d[strand][ref].keys():
                                
                                raw_coverage = p_db_raw.d[strand][ref][coord]
                                stat_dict['n_of_nts_covered']
                                
                                #construct the row tuple
                                t = (ref, coord, strand, raw_coverage, 0, 0)
                                cnt += 1
                                
                                # periodically report writeback tally
                                if cnt % 1000000 == 0:
                                    print('{0} writes to sqlite'.format(str(cnt)))
                                
                                yield t
                
                # inject sites into sqlite
                self.c.executemany('''INSERT INTO clip_signature VALUES
                                       (?,?,?,?,?,?)''', row_gen())
                
        
        else:
            # since no cleav_db_fp was provided, it is assumed that all
            # all reads are adapter-clipped
            self.c.execute('''ATTACH ? AS sam''',(sam_db_fp,))
            #select cleaved reads
            if uniq_only:
                # to select reads with non-null uniq_flag here
                self.c.execute('''SELECT ref TEXT, 
                                         start_coord INT, 
                                         end_coord INT, 
                                         strand TEXT, 
                                         cigar TEXT 
                                  FROM sam
                                  WHERE AND sam.uniq_flag NOT NULL''')
            else:
                self.c.execute('''SELECT ref TEXT, 
                                         start_coord INT, 
                                         end_coord INT, 
                                         strand TEXT, 
                                         cigar TEXT 
                                  FROM sam''')
            
            print('Tallying adapter-clipped reads')
            with PositionalShelf(os.path.join(temp_dir.name, sam_db_fn_prefix) + '.raw') as p_db_raw, \
                 PositionalShelf(os.path.join(temp_dir.name, sam_db_fn_prefix) + '.cleav') as p_db_cleav, \
                 PositionalShelf(os.path.join(temp_dir.name, sam_db_fn_prefix) + '.oneD') as p_db_oneD:
                
                n_of_alignments = 0
                l = self.c.fetchmany()
                while l:
                    for row in l:
                        stat_dict['n_of_mapped_cleav_reads'] += 1
                        if stat_dict['n_of_mapped_cleav_reads'] % 1000000 == 0:
                            print(stat_dict['n_of_mapped_cleav_reads'], '''alignments of adapter-clipped reads have been processed''')
                        
                        ref, start_coord, end_coord, strand, cigar = row
                        p_db_raw.addSiteTup( (ref, start_coord, end_coord, strand, 1) )
                        
                        #add left cleavage site
                        p_db_cleav.addSiteTup( (ref, start_coord, start_coord + 1, strand, 1) )
                        #add right cleavage site
                        p_db_cleav.addSiteTup( (ref, end_coord - 1, end_coord, strand, 1) )
                        
                        if '1D' in cigar:
                            oneDsite_l = get1DsiteList(ref, start_coord, strand, cigar)
                            if oneDsite_l:
                                stat_dict['n_of_mapped_oneD_reads'] += 1
                                for oneDsite in oneDsite_l:
                                    ref, start_coord, end_coord, strand = oneDsite
                                    p_db_oneD.addSiteTup( (ref, start_coord, end_coord, strand, 1) )
                    
                    l = self.c.fetchmany()
                
                print(stat_dict['n_of_mapped_cleav_reads'], '''alignments of adapter-clipped reads have been processed''')
                
                self.c.execute('''DETACH sam''')
                #TODO: remove cleav_db_fp...or put it in a self-destructing temp folder
                
                #perform data sync before write-out
                p_db_raw.d.sync()
                p_db_cleav.d.sync()
                p_db_oneD.d.sync()
                
                #TODO: optimize order of PRIMARY KEY indices
                self.c.execute('''CREATE TABLE clip_signature
                                     (
                                        ref TEXT, 
                                        coord INT, 
                                        strand TEXT, 
                                        raw_coverage INT, 
                                        cleavage INT, 
                                        single_nt_del INT, 
                                        PRIMARY KEY(strand, ref, coord)
                                     )''')
                
                # iterate over 1D sites first, lookup sites in cleav and raw
                # positional db and delete once they are written to the hits-clip 
                # sqlite to maintain unique site index across the shelves.
                def row_gen():
                    print('injecting clip signature values for 1D sites')
                    cnt = 0
                    for strand in p_db_oneD.d.keys():
                        print('loading sites from the {0} strand'.format(strand))
                        for ref in p_db_oneD.d[strand].keys():
                            for coord in p_db_oneD.d[strand][ref].keys():
                                
                                oneD_coverage = p_db_oneD.d[strand][ref][coord]
                                stat_dict['n_of_oneD_sites'] += 1
                                
                                # pull raw coverage value
                                raw_coverage = p_db_raw.d[strand][ref][coord]
                                del p_db_raw.d[strand][ref][coord]
                                stat_dict['n_of_nts_covered'] += 1
                                
                                # pull cleav coverage value
                                if p_db_cleav.key_exists( ref, strand, coord):
                                    cleav_coverage = p_db_cleav.d[strand][ref][coord]
                                    del p_db_cleav.d[strand][ref][coord]
                                    stat_dict['n_of_cleavage_sites'] += 1
                                else:
                                    cleav_coverage = 0
                                
                                #construct the row tuple
                                t = (ref, coord, strand, raw_coverage, cleav_coverage, oneD_coverage)
                                cnt += 1
                                
                                # periodically sync data to writeback key removals
                                if cnt % 1000000 == 0:
                                    p_db_raw.d.sync()
                                    p_db_cleav.d.sync()
                                    print('{0} writes to sqlite'.format(str(cnt)))
                                
                                yield t
                
                # inject sites into sqlite
                self.c.executemany('''INSERT INTO clip_signature VALUES
                                       (?,?,?,?,?,?)''', row_gen())
                
                
                # sync positional shelf data to writeback key removals
                p_db_raw.d.sync()
                p_db_cleav.d.sync()
                
                # iterate over cleaved sites, lookup sites in raw
                # positional db and delete once they are written to the hits-clip 
                # sqlite to maintain unique site index across the shelves.
                def row_gen():
                    print('injecting clip signature values for remaining cleavage sites')
                    cnt = 0
                    for strand in p_db_cleav.d.keys():
                        print('loading sites from the {0} strand'.format(strand))
                        for ref in p_db_cleav.d[strand].keys():
                            for coord in p_db_cleav.d[strand][ref].keys():
                                
                                cleav_coverage = p_db_cleav.d[strand][ref][coord]
                                stat_dict['n_of_cleavage_sites'] += 1
                                
                                # pull raw coverage value
                                raw_coverage = p_db_raw.d[strand][ref][coord]
                                del p_db_raw.d[strand][ref][coord]
                                stat_dict['n_of_nts_covered']
                                
                                #construct the row tuple
                                t = (ref, coord, strand, raw_coverage, cleav_coverage, 0)
                                cnt += 1
                                
                                # periodically sync data to writeback key removals
                                if cnt % 1000000 == 0:
                                    p_db_raw.d.sync()
                                    print('{0} writes to sqlite'.format(str(cnt)))
                                
                                yield t
                
                # inject sites into sqlite
                self.c.executemany('''INSERT INTO clip_signature VALUES
                                       (?,?,?,?,?,?)''', row_gen())
                
                
                # sync positional shelf data to writeback key removals
                p_db_raw.d.sync()
                
                # iterate over the remaining sites
                def row_gen():
                    print('injecting remaining raw coverage sites')
                    cnt = 0
                    for strand in p_db_raw.d.keys():
                        print('loading sites from the {0} strand'.format(strand))
                        for ref in p_db_raw.d[strand].keys():
                            for coord in p_db_raw.d[strand][ref].keys():
                                
                                raw_coverage = p_db_raw.d[strand][ref][coord]
                                stat_dict['n_of_nts_covered'] += 1
                                
                                #construct the row tuple
                                t = (ref, coord, strand, raw_coverage, 0, 0)
                                cnt += 1
                                
                                # periodically report writeback tally
                                if cnt % 1000000 == 0:
                                    print('{0} writes to sqlite'.format(str(cnt)))
                                
                                yield t
                
                # inject sites into sqlite
                self.c.executemany('''INSERT INTO clip_signature VALUES
                                       (?,?,?,?,?,?)''', row_gen())
                
        
        # Remove positional shelf files, THEY ARE WORTHLESS because entries
        # are deleted along the way.
        temp_dir.cleanup()
        
        self.conn.commit()
        
        return stat_dict
    
    def siteTup2sigTup(self, 
                       siteTup, 
                       flank = 0, 
                       cleav_rateMode = False, 
                       oneD_rateMode = True, 
                       rateMode_cut = 15, 
                       stranded = True):
        #HITSCLIP and iCLIP via inheritance
        """
        A method to query the db with a siteTup. This method will
        return a tuple of 3 lists, that give basewise coverage values 
        encompassed by the coords supplied in siteTup.
        siteTup = ( reference, startCoord, endCoord, strand )
        NOTE: '-' strand data is reversed so that all sigTups
        contain arrays that are indexed from 5' to 3'
        """
        #TODO: If the data from both strands should be aggregated 
        #      (as in strand-unaware seq technologies) use 'both'
        #      as the strand identifier.) NOTE: NOT IMPLEMENTED YET
        
        q_ref, q_startCoord, q_endCoord, q_strand = siteTup
        q_startCoord -= flank
        q_endCoord += flank
        t = (q_ref, q_strand, q_startCoord, q_endCoord)
        
        if q_strand != '+' and q_strand != '-':
            print('Cannot use strand info for strand type {0}'.format(strand))
            raise IOError
            #TODO: support datasets with no strand info (e.g. rip-seq) using "both"
            #this involves simply adding together the + and - strand entries
        
        self.c.execute('''SELECT * FROM clip_signature WHERE ref = ? AND strand = ? AND coord BETWEEN ?-1 AND ?''', t)
        dataD = {}
        for db_ref, db_coord, db_strand, raw_coverage, cleavage, single_nt_del in self.c.fetchall():
            dataD[db_coord] = (raw_coverage, cleavage, single_nt_del)
        
        rawCovL = []
        cleavCovL = []
        oneDcovL = []
        
        
        for i in range(q_startCoord, q_endCoord):
            rawCov, cleavCov, oneDcov = dataD.get(i, (0,0,0))
            rawCovL.append(rawCov)
            cleavCovL.append(cleavCov)
            oneDcovL.append(oneDcov)
        
        if not stranded:
            if q_strand == '+': 
                other_strand = '-'
            elif q_strand == '-': 
                other_strand = '+'
            
            t = (q_ref, other_strand, q_startCoord, q_endCoord)
            self.c.execute('''SELECT * FROM clip_signature WHERE ref = ? AND strand = ? AND coord BETWEEN ?-1 AND ?''', t)
            dataD = {}
            for db_ref, db_coord, db_strand, raw_coverage, cleavage, single_nt_del in self.c.fetchall():
                dataD[db_coord] = (raw_coverage, cleavage, single_nt_del)
            
            for idx, i in enumerate(range(q_startCoord, q_endCoord)):
                rawCov, cleavCov, oneDcov = dataD.get(i, (0,0,0))
                rawCovL[idx] += rawCov
                cleavCovL[idx] += cleavCov
                oneDcovL[idx] += oneDcov
        
        
        #put everything in the same orientation
        if q_strand == '-':
            rawCovL.reverse()
            cleavCovL.reverse()
            oneDcovL.reverse()
        
        if oneD_rateMode == True:
            oneDrateL = []
            for x, y in zip(oneDcovL, rawCovL):
                if y < rateMode_cut:
                    oneDrateL.append(0)
                else:
                    oneDrateL.append(x/y)
            
            oneDcovL = oneDrateL
        
        if cleav_rateMode == True:
            cleavRateL = []
            for x, y in zip(cleavCovL, rawCovL):
                if y < rateMode_cut:
                    cleavRateL.append(0)
                else:
                    r = x/y
                    cleavRateL.append(x/y)
            
            cleavCovL = cleavRateL
        
        return (rawCovL, cleavCovL, oneDcovL)

#TODO: incorporate the oneD (or other CIMS) site permutation routine
##### if you decide to store the oneD data you can use struct to generate byte
##### arrays encoding the one d counts for each row and give a NULL if no
##### oneD is present

class IclipSQLite(HitsClipSQLite):
    """
    SQLite adapters for iClip data.
    """
    
    def build_basewise_signature(self, sam_db_fp, uniq_only = False):
        '''
        Generate basewise counts for 1D, cleavage and raw coverage from the 
        set of alignements at the bwa-samse output file at sam_db_fp.
        The build process first parses the sam file and exports the calculations
        to positional shelve db files at pos_shelve_fp. The export to shelf 
        is necessary because building the positional sqlite database takes too 
        much time due to the need for UPDATE calls to read each new site. 
        Therefore, the PositonalShelves are written and then copied 
        to the connected sqlite database.
        
        Note: if the cleav_db_fp is set to None then it is assumed that all
        the reads in the sam file are adapater-clipped
        '''
        
        #NOTE: THE PRIMARY DIFFERENCE IS THAT ALL THE 5' TERMINI ARE
        # COUNTED IN THE CLEAVAGE COVERAGE VECTOR.
        
        #TODO: sqlite file with sam and fq tables can be temporary
        
        import os
        import tempfile
        import re
        
        sam_db_fn_prefix = os.path.basename(sam_db_fp).split('.')[0]
        temp_dir = tempfile.TemporaryDirectory( suffix=sam_db_fn_prefix, 
                                                dir=os.path.dirname(self.fp) )
        
        CIGARregexP = re.compile(r'([0-9]|[0-9][0-9])([MIDNSHP=X])')
        
        def get1DsiteList(ref, start_coord, strand, cigar_string):
            '''
            A function to glean the list of site tuples comprised of all
            1D events described by the zero-based start_coord argument and
            the cigar string that was copied from the sam file and is
            therefore one-based.
            '''
            oneDSitesList = []
            m = CIGARregexP.findall(cigar_string)
            for idx, tup in enumerate(m):
                posCount, operation = tup
                if operation == 'D':
                    offset = sum([int(t[0]) for t in m[:idx + 1] if t[1] != 'I'])
                    offset -= 1 # convert to zero-based coord
                    oneDsiteTup = (ref, 
                                   start_coord + offset, 
                                   start_coord + offset + 1, 
                                   strand)
                    oneDSitesList.append(oneDsiteTup)
            return oneDSitesList
        
        stat_dict = {}
        stat_dict['n_of_nts_covered'] = 0 # total of nt sites with coverage
        stat_dict['n_of_oneD_sites'] = 0 # total number of nts with oneD
        stat_dict['n_of_mapped_oneD_reads'] = 0 #note: impacted by uniq_only
        stat_dict['n_of_reads_processed'] = 0
        stat_dict['n_of_cleavage_sites'] = 0
        
        # since no cleav_db_fp was provided, it is assumed that all
        # all reads are adapter-clipped
        self.c.execute('''ATTACH ? AS sam''',(sam_db_fp,))
        #select cleaved reads
        if uniq_only:
            # to select reads with non-null uniq_flag here
            self.c.execute('''SELECT ref TEXT, 
                                     start_coord INT, 
                                     end_coord INT, 
                                     strand TEXT, 
                                     cigar TEXT 
                              FROM sam
                              WHERE sam.uniq_flag NOT NULL''')
        else:
            self.c.execute('''SELECT ref TEXT, 
                                     start_coord INT, 
                                     end_coord INT, 
                                     strand TEXT, 
                                     cigar TEXT 
                              FROM sam''')
        
        with PositionalShelf(os.path.join(temp_dir.name, sam_db_fn_prefix) + '.raw') as p_db_raw, \
             PositionalShelf(os.path.join(temp_dir.name, sam_db_fn_prefix) + '.cleav') as p_db_cleav, \
             PositionalShelf(os.path.join(temp_dir.name, sam_db_fn_prefix) + '.oneD') as p_db_oneD:
            
            l = self.c.fetchmany()
            while l:
                for row in l:
                    stat_dict['n_of_reads_processed'] += 1
                    if stat_dict['n_of_reads_processed'] % 1000000 == 0:
                        print(stat_dict['n_of_reads_processed'], '''alignments have been processed''')
                    
                    ref, start_coord, end_coord, strand, cigar = row
                    p_db_raw.addSiteTup( (ref, start_coord, end_coord, strand, 1) )
                    
                    #add left cleavage site
                    p_db_cleav.addSiteTup( (ref, start_coord, start_coord + 1, strand, 1) )
                    
                    if '1D' in cigar:
                        oneDsite_l = get1DsiteList(ref, start_coord, strand, cigar)
                        if oneDsite_l:
                            stat_dict['n_of_mapped_oneD_reads'] += 1
                            for oneDsite in oneDsite_l:
                                ref, start_coord, end_coord, strand = oneDsite
                                p_db_oneD.addSiteTup( (ref, start_coord, end_coord, strand, 1) )
                
                l = self.c.fetchmany()
            
            print(stat_dict['n_of_reads_processed'], '''alignments have been processed''')
            
            self.c.execute('''DETACH sam''')
            
            #perform data sync before write-out
            p_db_raw.d.sync()
            p_db_cleav.d.sync()
            p_db_oneD.d.sync()
            
            #TODO: optimize order of PRIMARY KEY indices
            self.c.execute('''CREATE TABLE clip_signature
                                 (
                                    ref TEXT, 
                                    coord INT, 
                                    strand TEXT, 
                                    raw_coverage INT, 
                                    cleavage INT, 
                                    single_nt_del INT, 
                                    PRIMARY KEY(strand, ref, coord)
                                 )''')
            
            # iterate over 1D sites first, lookup sites in cleav and raw
            # positional db and delete once they are written to the hits-clip 
            # sqlite to maintain unique site index across the shelves.
            def row_gen():
                print('injecting clip signature values for 1D sites')
                cnt = 0
                for strand in p_db_oneD.d.keys():
                    print('loading sites from the {0} strand'.format(strand))
                    for ref in p_db_oneD.d[strand].keys():
                        for coord in p_db_oneD.d[strand][ref].keys():
                            
                            oneD_coverage = p_db_oneD.d[strand][ref][coord]
                            stat_dict['n_of_oneD_sites'] += 1
                            
                            # pull raw coverage value
                            raw_coverage = p_db_raw.d[strand][ref][coord]
                            del p_db_raw.d[strand][ref][coord]
                            stat_dict['n_of_nts_covered'] += 1
                            
                            # pull cleav coverage value
                            if p_db_cleav.key_exists( ref, strand, coord):
                                cleav_coverage = p_db_cleav.d[strand][ref][coord]
                                del p_db_cleav.d[strand][ref][coord]
                                stat_dict['n_of_cleavage_sites'] += 1
                            else:
                                cleav_coverage = 0
                            
                            #construct the row tuple
                            t = (ref, coord, strand, raw_coverage, cleav_coverage, oneD_coverage)
                            cnt += 1
                            
                            # periodically sync data to writeback key removals
                            if cnt % 1000000 == 0:
                                p_db_raw.d.sync()
                                p_db_cleav.d.sync()
                                print('{0} writes to sqlite'.format(str(cnt)))
                            
                            yield t
            
            # inject sites into sqlite
            self.c.executemany('''INSERT INTO clip_signature VALUES
                                   (?,?,?,?,?,?)''', row_gen())
            
            
            # sync positional shelf data to writeback key removals
            p_db_raw.d.sync()
            p_db_cleav.d.sync()
            
            # iterate over cleaved sites, lookup sites in raw
            # positional db and delete once they are written to the hits-clip 
            # sqlite to maintain unique site index across the shelves.
            def row_gen():
                print('injecting clip signature values for remaining cleavage sites')
                cnt = 0
                for strand in p_db_cleav.d.keys():
                    print('loading sites from the {0} strand'.format(strand))
                    for ref in p_db_cleav.d[strand].keys():
                        for coord in p_db_cleav.d[strand][ref].keys():
                            
                            cleav_coverage = p_db_cleav.d[strand][ref][coord]
                            stat_dict['n_of_cleavage_sites'] += 1
                            
                            # pull raw coverage value
                            raw_coverage = p_db_raw.d[strand][ref][coord]
                            del p_db_raw.d[strand][ref][coord]
                            stat_dict['n_of_nts_covered']
                            
                            #construct the row tuple
                            t = (ref, coord, strand, raw_coverage, cleav_coverage, 0)
                            cnt += 1
                            
                            # periodically sync data to writeback key removals
                            if cnt % 1000000 == 0:
                                p_db_raw.d.sync()
                                print('{0} writes to sqlite'.format(str(cnt)))
                            
                            yield t
            
            # inject sites into sqlite
            self.c.executemany('''INSERT INTO clip_signature VALUES
                                   (?,?,?,?,?,?)''', row_gen())
            
            
            # sync positional shelf data to writeback key removals
            p_db_raw.d.sync()
            
            # iterate over the remaining sites
            def row_gen():
                print('injecting remaining raw coverage sites')
                cnt = 0
                for strand in p_db_raw.d.keys():
                    print('loading sites from the {0} strand'.format(strand))
                    for ref in p_db_raw.d[strand].keys():
                        for coord in p_db_raw.d[strand][ref].keys():
                            
                            raw_coverage = p_db_raw.d[strand][ref][coord]
                            stat_dict['n_of_nts_covered'] += 1
                            
                            #construct the row tuple
                            t = (ref, coord, strand, raw_coverage, 0, 0)
                            cnt += 1
                            
                            # periodically report writeback tally
                            if cnt % 1000000 == 0:
                                print('{0} writes to sqlite'.format(str(cnt)))
                            
                            yield t
            
            # inject sites into sqlite
            self.c.executemany('''INSERT INTO clip_signature VALUES
                                   (?,?,?,?,?,?)''', row_gen())
        
        
        # Remove positional shelf files, HEY ARE WORTHLESS because entries
        # are deleted along the way.
        temp_dir.cleanup()
        
        self.conn.commit()
        
        return stat_dict



class ParClipSQLite(HitsClipSQLite):
    """
    SQLite adapters for ParClip data.
    """
    
    def input_sam(self, fp, uniq_only = False):
        """
        Input a sam file.
        """
        #NOTE: If you filter for uniq here you have a much smaller db and you save time.
        #HOWEVER! if you are going to do a permutation routine downstream then 
        #      you want to set uniq_only to false and use the uniq_only selection
        #      downstream (in build_basewise_signature) so that all the 1D events 
        #      are available in the sam sqlite for permutations
        
        #TODO: include the function to parse MAXOCC flag in bwa alignments
        #      that were invoked with the maxOcc option defined
        from clipPyL.flatfile_parsing import SamReader
        
        self.c.execute('''CREATE TABLE sam
                     (
                        qname TEXT, 
                        ref TEXT, 
                        start_coord INT, 
                        end_coord INT, 
                        strand TEXT, 
                        t_to_c_coords TEXT, 
                        uniq_flag INT,
                        PRIMARY KEY(qname)
                     )''')
        
        stat_dict = {}
        stat_dict['n_of_reads'] = 0 #total number of reads is gleaned from sam file
        stat_dict['n_of_mapped_reads'] = 0 #gleaned from sam flagbits
        stat_dict['n_of_unique_aligns'] = 0 #gleaned from sam XT opt field
        
        with SamReader(fp) as sam_gen:
            
            def row_gen():
                
                for d in sam_gen:
                    stat_dict['n_of_reads'] += 1
                    
                    if stat_dict['n_of_reads'] % 1000000 == 0:
                        print(stat_dict['n_of_reads'], ' reads processed')
                    
                    if sam_gen.isMapped():
                        stat_dict['n_of_mapped_reads'] += 1
                        
                        #identify uniquely aligned reads by accessing 
                        #the sam optdict
                        if sam_gen.unique_bool():
                            stat_dict['n_of_unique_aligns'] += 1
                            uniq_flag = 1
                        elif uniq_only:
                            continue
                        else:
                            uniq_flag = None
                        
                        #raw coverage interval
                        reference, startCoord, endCoord, strand = sam_gen.calc_ucsc_interval()
                        
                        # detect the t-to-c mutations and store the zero-based reference 
                        # coordinate of each in a comma-seperated string (will be parsed
                        # downstream in build_basewise_signature. othewise, NULL
                        if d['OPTDICT']:
                            if strand == '+':
                                v = sam_gen.detect_mm_variant(native_nt = 'T', variant_nt = 'C')
                            elif strand == '-':
                                #NOTE: this is necessary because reverse complement of T->C is A->G
                                v = sam_gen.detect_mm_variant(native_nt = 'A', variant_nt = 'G')
                            
                            if v:
                                t_to_c_coordL = []
                                for sam_seq_pos, ref_pos in v:
                                    t_to_c_coordL.append(str(ref_pos))
                                yield (d['QNAME'], reference, startCoord, endCoord, strand, ','.join(t_to_c_coordL) , uniq_flag)
                            else:
                                yield (d['QNAME'], reference, startCoord, endCoord, strand, None, uniq_flag)
                        else:
                            yield (d['QNAME'], reference, startCoord, endCoord, strand, None, uniq_flag)
            
            self.c.executemany('''INSERT INTO sam VALUES (?,?,?,?,?,?,?)''', row_gen())
        
        self.conn.commit()
        
        return stat_dict
    
    def build_basewise_signature(self, sam_db_fp, cleav_db_fp = None, uniq_only = False, chunk = 10000):
        '''
        Generate basewise counts for T-to-C, cleavage and raw coverage from the 
        set of alignements at the bwa-samse output file at sam_db_fp.
        The build process first parses the sam file and exports the calculations
        to positional shelve db files at pos_shelve_fp. The export to shelf 
        is necessary because building the positional sqlite database takes too 
        much time due to the need for UPDATE calls to read each new site. 
        Therefore, the PositonalShelves are written and then copied 
        to the connected sqlite database.
        
        Note: if the cleav_db_fp is set to None then it is assumed that all
        the reads in the sam file are adapater-clipped
        '''
        #TODO: positional db can be a temporary file
        #TODO: sqlite file with sam and fq tables can be temporary
        # the way to do this will be to use temp directory file paths
        # in context managers to generate fp arguments that are passed to 
        # these functions
        #TODO: add an option to input a simple list of clipped reads upstream 
        # at "input_fastq" method. ...which should be changed 
        # to "build_clipped_readid_db" method
        import os
        import tempfile
        import re
        
        sam_db_fn_prefix = os.path.basename(sam_db_fp).split('.')[0]
        temp_dir = tempfile.TemporaryDirectory( suffix=sam_db_fn_prefix, 
                                                dir=os.path.dirname(self.fp) )
        
        #LB_snip
        
        def getTtoCsiteList(ref, start_coord, strand, t_to_c_coordL): #cigar_string):
            '''
            A function to glean the list of site tuples comprised of all
            T->C events described by the zero-based start_coord argument and
            the t_to_c_coodL.
            '''
            t_to_c_sitesL = []
            l = t_to_c_coordL.split(',')
            for coord in l:
                if coord:
                    t_to_c_siteTup = (ref, 
                                      int(coord), 
                                      int(coord) + 1, 
                                      strand)
                    t_to_c_sitesL.append(t_to_c_siteTup)
            
            return t_to_c_sitesL
        
        stat_dict = {}
        stat_dict['n_of_nts_covered'] = 0 # total of nt sites with coverage
        stat_dict['n_of_t_to_c_sites'] = 0 # total number of nts with t to c substitution
        stat_dict['n_of_cleavage_sites'] = 0 #total number of nts with cleavage
        stat_dict['n_of_mapped_cleav_reads'] = 0 #note: impacted by uniq_only
        stat_dict['n_of_mapped_t_to_c_reads'] = 0 #note: impacted by uniq_only
        
        if cleav_db_fp:
            
            print('Selecting adapter-clipped reads')
            self.c.execute('''ATTACH ? AS adpt_clipped_reads''',(cleav_db_fp,))
            self.c.execute('''ATTACH ? AS sam''',(sam_db_fp,))
            #select cleaved reads
            if uniq_only:
                # to select reads with non-null uniq_flag here
                self.c.execute('''SELECT ref TEXT, 
                                         start_coord INT, 
                                         end_coord INT, 
                                         strand TEXT, 
                                         t_to_c_coords TEXT 
                                  FROM sam, adpt_clipped_reads 
                                  WHERE sam.qname = adpt_clipped_reads.read_id 
                                  AND sam.uniq_flag NOT NULL''')
            else:
                self.c.execute('''SELECT ref TEXT, 
                                         start_coord INT, 
                                         end_coord INT, 
                                         strand TEXT, 
                                         t_to_c_coords TEXT 
                                  FROM sam, adpt_clipped_reads 
                                  WHERE sam.qname = adpt_clipped_reads.read_id''')
            
            print('Tallying adapter-clipped reads')
            # subdirectory and explicity delete it in the script
            with PositionalShelf(os.path.join(temp_dir.name, sam_db_fn_prefix) + '.raw') as p_db_raw, \
                 PositionalShelf(os.path.join(temp_dir.name, sam_db_fn_prefix) + '.cleav') as p_db_cleav, \
                 PositionalShelf(os.path.join(temp_dir.name, sam_db_fn_prefix) + '.t_to_c') as p_db_t_to_c:
                
                l = self.c.fetchmany()
                while l:
                    for row in l:
                        stat_dict['n_of_mapped_cleav_reads'] += 1
                        if stat_dict['n_of_mapped_cleav_reads'] % 1000000 == 0:
                            print(stat_dict['n_of_mapped_cleav_reads'], '''alignments of adapter-clipped reads have been processed''')
                        
                        ref, start_coord, end_coord, strand, t_to_c_coords = row
                        p_db_raw.addSiteTup( (ref, start_coord, end_coord, strand, 1) )
                        
                        #add left cleavage site
                        p_db_cleav.addSiteTup( (ref, start_coord, start_coord + 1, strand, 1) )
                        #add right cleavage site
                        p_db_cleav.addSiteTup( (ref, end_coord - 1, end_coord, strand, 1) )
                        
                        if t_to_c_coords:
                            t_to_c_siteL = getTtoCsiteList(ref, start_coord, strand, t_to_c_coords)
                            if t_to_c_siteL:
                                stat_dict['n_of_mapped_t_to_c_reads'] += 1
                                for t_to_c_site in t_to_c_siteL:
                                    ref, start_coord, end_coord, strand = t_to_c_site
                                    p_db_t_to_c.addSiteTup( (ref, start_coord, end_coord, strand, 1) )
                    
                    l = self.c.fetchmany()
                
                #perform data sync before write-out
                p_db_raw.d.sync()
                p_db_cleav.d.sync()
                p_db_t_to_c.d.sync()
                
                print(stat_dict['n_of_mapped_cleav_reads'], '''alignments of adapter-clipped reads have been processed''')
                
                print('Selecting unclipped reads')
                if uniq_only:
                    self.c.execute('''SELECT ref TEXT, 
                                             start_coord INT, 
                                             end_coord INT, 
                                             strand TEXT, 
                                             t_to_c_coords TEXT 
                                      FROM sam
                                      LEFT JOIN adpt_clipped_reads ON adpt_clipped_reads.read_id = sam.qname
                                      WHERE adpt_clipped_reads.read_id IS NULL
                                      AND sam.uniq_flag NOT NULL''')
                else:
                    self.c.execute('''SELECT ref TEXT, 
                                             start_coord INT, 
                                             end_coord INT, 
                                             strand TEXT, 
                                             t_to_c_coords TEXT 
                                      FROM sam
                                      LEFT JOIN adpt_clipped_reads ON adpt_clipped_reads.read_id = sam.qname
                                      WHERE adpt_clipped_reads.read_id IS NULL''')
                
                print('Tallying unclipped reads.')
                n_of_alignments = 0
                l = self.c.fetchmany()
                while l:
                    for row in l:
                        n_of_alignments += 1
                        if n_of_alignments % 1000000 == 0:
                            print(n_of_alignments, '''alignments of unclipped reads have been processed''')
                        
                        ref, start_coord, end_coord, strand, t_to_c_coords = row
                        p_db_raw.addSiteTup( (ref, start_coord, end_coord, strand, 1) )
                        
                        if t_to_c_coords:
                            t_to_c_siteL = getTtoCsiteList(ref, start_coord, strand, t_to_c_coords)
                            if t_to_c_siteL:
                                stat_dict['n_of_mapped_t_to_c_reads'] += 1
                                for t_to_c_site in t_to_c_siteL:
                                    ref, start_coord, end_coord, strand = t_to_c_site
                                    p_db_t_to_c.addSiteTup( (ref, start_coord, end_coord, strand, 1) )
                    
                    l = self.c.fetchmany()
                
                print(n_of_alignments, '''alignments of unclipped reads have been processed''')
                
                self.c.execute('''DETACH adpt_clipped_reads''')
                self.c.execute('''DETACH sam''')
                
                #perform data sync before write-out
                p_db_raw.d.sync()
                p_db_cleav.d.sync()
                p_db_t_to_c.d.sync()
                
                #TODO: optimize order of PRIMARY KEY indices
                self.c.execute('''CREATE TABLE clip_signature
                                     (
                                        ref TEXT, 
                                        coord INT, 
                                        strand TEXT, 
                                        raw_coverage INT, 
                                        cleavage INT, 
                                        t_to_c INT, 
                                        PRIMARY KEY(strand, ref, coord)
                                     )''')
                
                # iterate over t_to_c sites first, lookup sites in cleav and raw
                # positional db and delete once they are written to the hits-clip 
                # sqlite to maintain unique site index across the shelves.
                def row_gen():
                    print('injecting clip signature values for t_to_c sites')
                    cnt = 0
                    for strand in p_db_t_to_c.d.keys():
                        print('loading sites from the {0} strand'.format(strand))
                        for ref in p_db_t_to_c.d[strand].keys():
                            for coord in p_db_t_to_c.d[strand][ref].keys():
                                
                                t_to_c_coverage = p_db_t_to_c.d[strand][ref][coord]
                                stat_dict['n_of_t_to_c_sites'] += 1
                                
                                # pull raw coverage value
                                raw_coverage = p_db_raw.d[strand][ref][coord]
                                del p_db_raw.d[strand][ref][coord]
                                stat_dict['n_of_nts_covered'] += 1
                                
                                # pull cleav coverage value
                                if p_db_cleav.key_exists( ref, strand, coord):
                                    cleav_coverage = p_db_cleav.d[strand][ref][coord]
                                    del p_db_cleav.d[strand][ref][coord]
                                    stat_dict['n_of_cleavage_sites'] += 1
                                else:
                                    cleav_coverage = 0
                                
                                #construct the row tuple
                                t = (ref, coord, strand, raw_coverage, cleav_coverage, t_to_c_coverage)
                                cnt += 1
                                
                                # periodically sync data to writeback key removals
                                if cnt % chunk == 0:
                                    p_db_raw.d.sync()
                                    p_db_cleav.d.sync()
                                    print('{0} writes to sqlite'.format(str(cnt)))
                                
                                yield t
                
                # inject sites into sqlite
                self.c.executemany('''INSERT INTO clip_signature VALUES
                                       (?,?,?,?,?,?)''', row_gen())
                
                # sync positional shelf data to writeback key removals
                p_db_raw.d.sync()
                p_db_cleav.d.sync()
                
                # iterate over the remaining cleaved sites, lookup sites in raw
                # positional db and delete once they are written to the hits-clip 
                # sqlite to maintain unique site index across the shelves.
                def row_gen():
                    print('injecting clip signature values for remaining cleavage sites')
                    cnt = 0
                    for strand in p_db_cleav.d.keys():
                        print('loading sites from the {0} strand'.format(strand))
                        for ref in p_db_cleav.d[strand].keys():
                            for coord in p_db_cleav.d[strand][ref].keys():
                                
                                cleav_coverage = p_db_cleav.d[strand][ref][coord]
                                stat_dict['n_of_cleavage_sites'] += 1
                                
                                # pull raw coverage value
                                raw_coverage = p_db_raw.d[strand][ref][coord]
                                del p_db_raw.d[strand][ref][coord]
                                stat_dict['n_of_nts_covered'] += 1
                                
                                #construct the row tuple
                                t = (ref, coord, strand, raw_coverage, cleav_coverage, 0)
                                cnt += 1
                                
                                # periodically sync data to writeback key removals
                                if cnt % chunk == 0:
                                    p_db_raw.d.sync()
                                    print('{0} writes to sqlite'.format(str(cnt)))
                                
                                yield t
                
                # inject sites into sqlite
                self.c.executemany('''INSERT INTO clip_signature VALUES
                                       (?,?,?,?,?,?)''', row_gen())
                
                # sync positional shelf data to writeback key removals
                p_db_raw.d.sync()
                
                # iterate over the remaining sites
                def row_gen():
                    print('injecting raw coverage values for remaining sites')
                    cnt = 0
                    for strand in p_db_raw.d.keys():
                        print('loading sites from the {0} strand'.format(strand))
                        for ref in p_db_raw.d[strand].keys():
                            for coord in p_db_raw.d[strand][ref].keys():
                                
                                raw_coverage = p_db_raw.d[strand][ref][coord]
                                stat_dict['n_of_nts_covered']
                                
                                #construct the row tuple
                                t = (ref, coord, strand, raw_coverage, 0, 0)
                                cnt += 1
                                
                                # periodically report writeback tally
                                if cnt % chunk == 0:
                                    print('{0} writes to sqlite'.format(str(cnt)))
                                
                                yield t
                
                # inject sites into sqlite
                self.c.executemany('''INSERT INTO clip_signature VALUES
                                       (?,?,?,?,?,?)''', row_gen())
                
        
        else:
            # since no cleav_db_fp was provided, it is assumed that all
            # all reads are adapter-clipped
            self.c.execute('''ATTACH ? AS sam''',(sam_db_fp,))
            #select cleaved reads
            if uniq_only:
                # to select reads with non-null uniq_flag here
                self.c.execute('''SELECT ref TEXT, 
                                         start_coord INT, 
                                         end_coord INT, 
                                         strand TEXT, 
                                         t_to_c_coords TEXT 
                                  FROM sam, adpt_clipped_reads 
                                  WHERE sam.qname = adpt_clipped_reads.read_id 
                                  AND sam.uniq_flag NOT NULL''')
            else:
                self.c.execute('''SELECT ref TEXT, 
                                         start_coord INT, 
                                         end_coord INT, 
                                         strand TEXT, 
                                         t_to_c_coords TEXT 
                                  FROM sam, adpt_clipped_reads 
                                  WHERE sam.qname = adpt_clipped_reads.read_id''')
            
            print('Tallying adapter-clipped reads')
            # subdirectory and explicity delete it in the script
            with PositionalShelf(os.path.join(temp_dir.name, sam_db_fn_prefix) + '.raw') as p_db_raw, \
                 PositionalShelf(os.path.join(temp_dir.name, sam_db_fn_prefix) + '.cleav') as p_db_cleav, \
                 PositionalShelf(os.path.join(temp_dir.name, sam_db_fn_prefix) + '.t_to_c') as p_db_t_to_c:
                
                l = self.c.fetchmany()
                while l:
                    for row in l:
                        stat_dict['n_of_mapped_cleav_reads'] += 1
                        if stat_dict['n_of_mapped_cleav_reads'] % 1000000 == 0:
                            print(stat_dict['n_of_mapped_cleav_reads'], '''alignments of adapter-clipped reads have been processed''')
                        
                        ref, start_coord, end_coord, strand, t_to_c_coords = row
                        p_db_raw.addSiteTup( (ref, start_coord, end_coord, strand, 1) )
                        
                        #add left cleavage site
                        p_db_cleav.addSiteTup( (ref, start_coord, start_coord + 1, strand, 1) )
                        #add right cleavage site
                        p_db_cleav.addSiteTup( (ref, end_coord - 1, end_coord, strand, 1) )
                        
                        if t_to_c_coords:
                            t_to_c_siteL = getTtoCsiteList(ref, start_coord, strand, t_to_c_coords)
                            if t_to_c_siteL:
                                stat_dict['n_of_mapped_t_to_c_reads'] += 1
                                for t_to_c_site in t_to_c_siteL:
                                    ref, start_coord, end_coord, strand = t_to_c_site
                                    p_db_t_to_c.addSiteTup( (ref, start_coord, end_coord, strand, 1) )
                    
                    l = self.c.fetchmany()
                
                print(stat_dict['n_of_mapped_cleav_reads'], '''alignments of adapter-clipped reads have been processed''')
                
                self.c.execute('''DETACH sam''')
                #TODO: remove cleav_db_fp...or put it in a self-destructing temp folder
                
                #perform data sync before write-out
                p_db_raw.d.sync()
                p_db_cleav.d.sync()
                p_db_t_to_c.d.sync()
                
                #TODO: optimize order of PRIMARY KEY indices
                self.c.execute('''CREATE TABLE clip_signature
                                     (
                                        ref TEXT, 
                                        coord INT, 
                                        strand TEXT, 
                                        raw_coverage INT, 
                                        cleavage INT, 
                                        t_to_c INT, 
                                        PRIMARY KEY(strand, ref, coord)
                                     )''')
                
                # iterate over t_to_c sites first, lookup sites in cleav and raw
                # positional db and delete once they are written to the hits-clip 
                # sqlite to maintain unique site index across the shelves.
                def row_gen():
                    print('injecting clip signature values for t_to_c sites')
                    cnt = 0
                    for strand in p_db_t_to_c.d.keys():
                        print('loading sites from the {0} strand'.format(strand))
                        for ref in p_db_t_to_c.d[strand].keys():
                            for coord in p_db_t_to_c.d[strand][ref].keys():
                                
                                t_to_c_coverage = p_db_t_to_c.d[strand][ref][coord]
                                stat_dict['n_of_t_to_c_sites'] += 1
                                
                                # pull raw coverage value
                                raw_coverage = p_db_raw.d[strand][ref][coord]
                                del p_db_raw.d[strand][ref][coord]
                                stat_dict['n_of_nts_covered'] += 1
                                
                                # pull cleav coverage value
                                if p_db_cleav.key_exists( ref, strand, coord):
                                    cleav_coverage = p_db_cleav.d[strand][ref][coord]
                                    del p_db_cleav.d[strand][ref][coord]
                                    stat_dict['n_of_cleavage_sites'] += 1
                                else:
                                    cleav_coverage = 0
                                
                                #construct the row tuple
                                t = (ref, coord, strand, raw_coverage, cleav_coverage, t_to_c_coverage)
                                cnt += 1
                                
                                # periodically sync data to writeback key removals
                                if cnt % chunk == 0:
                                    p_db_raw.d.sync()
                                    p_db_cleav.d.sync()
                                    print('{0} writes to sqlite'.format(str(cnt)))
                                
                                yield t
                
                # inject sites into sqlite
                self.c.executemany('''INSERT INTO clip_signature VALUES
                                       (?,?,?,?,?,?)''', row_gen())
                
                # sync positional shelf data to writeback key removals
                p_db_raw.d.sync()
                p_db_cleav.d.sync()
                
                # iterate over the remaining cleaved sites, lookup sites in raw
                # positional db and delete once they are written to the hits-clip 
                # sqlite to maintain unique site index across the shelves.
                def row_gen():
                    print('injecting clip signature values for remaining cleavage sites')
                    cnt = 0
                    for strand in p_db_cleav.d.keys():
                        print('loading sites from the {0} strand'.format(strand))
                        for ref in p_db_cleav.d[strand].keys():
                            for coord in p_db_cleav.d[strand][ref].keys():
                                
                                cleav_coverage = p_db_cleav.d[strand][ref][coord]
                                stat_dict['n_of_cleavage_sites'] += 1
                                
                                # pull raw coverage value
                                raw_coverage = p_db_raw.d[strand][ref][coord]
                                del p_db_raw.d[strand][ref][coord]
                                stat_dict['n_of_nts_covered'] += 1
                                
                                #construct the row tuple
                                t = (ref, coord, strand, raw_coverage, cleav_coverage, 0)
                                cnt += 1
                                
                                # periodically sync data to writeback key removals
                                if cnt % chunk == 0:
                                    p_db_raw.d.sync()
                                    print('{0} writes to sqlite'.format(str(cnt)))
                                
                                yield t
                
                # inject sites into sqlite
                self.c.executemany('''INSERT INTO clip_signature VALUES
                                       (?,?,?,?,?,?)''', row_gen())
                
                # sync positional shelf data to writeback key removals
                p_db_raw.d.sync()
                
                # iterate over the remaining sites
                def row_gen():
                    print('injecting raw coverage values for remaining sites')
                    cnt = 0
                    for strand in p_db_raw.d.keys():
                        print('loading sites from the {0} strand'.format(strand))
                        for ref in p_db_raw.d[strand].keys():
                            for coord in p_db_raw.d[strand][ref].keys():
                                
                                raw_coverage = p_db_raw.d[strand][ref][coord]
                                stat_dict['n_of_nts_covered']
                                
                                #construct the row tuple
                                t = (ref, coord, strand, raw_coverage, 0, 0)
                                cnt += 1
                                
                                # periodically report writeback tally
                                if cnt % chunk == 0:
                                    print('{0} writes to sqlite'.format(str(cnt)))
                                
                                yield t
                
                # inject sites into sqlite
                self.c.executemany('''INSERT INTO clip_signature VALUES
                                       (?,?,?,?,?,?)''', row_gen())
                
        
        # Remove positional shelf files, HEY ARE WORTHLESS because entries
        # are deleted along the way.
        temp_dir.cleanup()
        
        self.conn.commit()
        
        return stat_dict
    
    def siteTup2sigTup(self, 
                       siteTup, 
                       flank = 0, 
                       t_to_c_rateMode = True, 
                       cleav_rateMode =False,
                       rateMode_cut = 15, 
                       stranded = True):
        #PARCLIP
        """
        A method to query the db with a siteTup. This method will
        return a tuple of 3 lists, that give basewise coverage values 
        encompassed by the coords supplied in siteTup.
        siteTup = ( reference, startCoord, endCoord, strand )
        NOTE: '-' strand data is reversed so that all sigTups
        contain arrays that are indexed from 5' to 3'
        """
        #TODO: If the data from both strands should be aggregated 
        #      (as in strand-unaware seq technologies) use 'both'
        #      as the strand identifier.) NOTE: NOT IMPLEMENTED YET
        
        q_ref, q_startCoord, q_endCoord, q_strand = siteTup
        q_startCoord -= flank
        q_endCoord += flank
        t = (q_ref, q_strand, q_startCoord, q_endCoord)
        
        if q_strand != '+' and q_strand != '-':
            print('Cannot use strand info for strand type {0}'.format(strand))
            raise IOError
            #TODO: support datasets with no strand info (e.g. rip-seq) using "both"
            #this involves simply adding together the + and - strand entries
        
        self.c.execute('''SELECT * FROM clip_signature WHERE ref = ? AND strand = ? AND coord BETWEEN ?-1 AND ?''', t)
        dataD = {}
        for db_ref, db_coord, db_strand, raw_coverage, cleavage, t_to_c in self.c.fetchall():
            dataD[db_coord] = (raw_coverage, cleavage, t_to_c)
        
        rawCovL = []
        cleavCovL = []
        t_to_c_covL = []
        
        
        for i in range(q_startCoord, q_endCoord):
            rawCov, cleavCov, t_to_c_cov = dataD.get(i, (0,0,0))
            rawCovL.append(rawCov)
            cleavCovL.append(cleavCov)
            t_to_c_covL.append(t_to_c_cov)
        
        if not stranded:
            if q_strand == '+': 
                other_strand = '-'
            elif q_strand == '-': 
                other_strand = '+'
            
            t = (q_ref, other_strand, q_startCoord, q_endCoord)
            self.c.execute('''SELECT * FROM clip_signature WHERE ref = ? AND strand = ? AND coord BETWEEN ?-1 AND ?''', t)
            dataD = {}
            for db_ref, db_coord, db_strand, raw_coverage, cleavage, t_to_c in self.c.fetchall():
                dataD[db_coord] = (raw_coverage, cleavage, t_to_c)
            
            for idx, i in enumerate(range(q_startCoord, q_endCoord)):
                rawCov, cleavCov, t_to_c_cov = dataD.get(i, (0,0,0))
                rawCovL[idx] += rawCov
                cleavCovL[idx] += cleavCov
                t_to_c_covL[idx] += t_to_c_cov
        
        
        #put everything in the same orientation
        if q_strand == '-':
            rawCovL.reverse()
            cleavCovL.reverse()
            t_to_c_covL.reverse()
        
        if t_to_c_rateMode == True:
            t_to_c_rateL = []
            for x, y in zip(t_to_c_covL, rawCovL):
                if y < rateMode_cut:
                    t_to_c_rateL.append(0)
                else:
                    t_to_c_rateL.append(x/y)
            
            t_to_c_covL = t_to_c_rateL
        
        if cleav_rateMode == True:
            cleavRateL = []
            for x, y in zip(cleavCovL, rawCovL):
                if y < rateMode_cut:
                    cleavRateL.append(0)
                else:
                    r = x/y
                    cleavRateL.append(x/y)
            
            cleavCovL = cleavRateL
        
        return (rawCovL, cleavCovL, t_to_c_covL)

class Bed6SQLite(SQLiteBase):
    
    def input_bed6(self, fp):
        
        self.c.execute('''CREATE TABLE bed6_intervals
                             (
                                chrom TEXT, 
                                chromStart INT, 
                                chromEnd INT, 
                                name TEXT, 
                                score INT, 
                                strand TEXT, 
                                PRIMARY KEY(strand, chrom, chromStart, chromEnd)
                             )''')
        
        from clipPyL.flatfile_parsing import Bed6Reader
        
        stat_dict = {}
        stat_dict['n_of_intervals'] = 0 #total number of intervals
        
        with Bed6Reader(fp) as bed6_gen:
            
            def row_gen():
                
                for d in bed6_gen:
                    
                    stat_dict['n_of_intervals'] += 1
                    
                    yield (d['chrom'], 
                           d['chromStart'], 
                           d['chromEnd'], 
                           d['name'], 
                           d['score'], 
                           d['strand'])
            
            self.c.executemany('''INSERT INTO bed6_intervals VALUES (?,?,?,?,?,?)''', row_gen())
        
        print('Read', stat_dict['n_of_intervals'], 'intervals')
        self.conn.commit()
        return
    
    def siteTup_lookup(self, query_siteTup):
        '''
        Look for all intervals that overlap with the query_siteTup.
        Return a list of dictionaries, one for each bed6 entry
        that intersects with the query_siteTup.
        '''
        # q_ prefix denotes "query"
        q_chrom, q_chromStart, q_chromEnd, q_strand = query_siteTup
        
        t = (q_strand, q_chrom, q_chromStart, q_chromEnd, q_chromStart, q_chromEnd)
        
        self.c.execute('''SELECT * FROM bed6_intervals
                            WHERE strand = ?
                            AND chrom = ?
                            AND (chromStart BETWEEN ? AND ? OR chromEnd BETWEEN ? AND ?)''', t)
        
        return self.c.fetchall()







