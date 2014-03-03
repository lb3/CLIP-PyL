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

#TO ACCOMODATE RNA ANNOTATIONS YOU SHOULD
#CREATE ANOTHER LEVEL OF BRANCHING IN THE
#SHELF SUCH THAT A RNA ANNOTATION ID CAN BE PARSED
#TO ITS CONSTITUENT PARTS TO SUPPLY KEYS TO NESTED LEVELS
#EX: ens098767800 can be parsed to ens key and then the integer can
# be parsed into ranges...this is likely necessary because the enormous amount
# of rna ids will cause slowdown during key lookup.

# an alternative approach is to just directly insert into the
# hitsclipsql db (see build_index_from_sam method)

#...or a acombination of the two where and additional
# text column is added and the index is built on
# id prefix and suffix columns such that 
# ens098767800 would become a text column and
# an integer column

class PositionalShelf():
    """
    API DEVELOPMENT WARNING: ALWAYS EXPLICITLY DECLARE BOTHE STRAND KEYS BECAUSE
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
    
    def __init__(self, fp, wb_bool = True, dataType = int, write_chunk = 100000):
        
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
        #check that the dataObject is expected type
        if type(dataObject) != self.dataType:
            print('incorrect data type given in dataObject component of the siteTup')
            raise TypeError
        
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
        elif self.dataType == str:
            #check that the site has been written to the chromDict in self.d storage object
            if reference in self.d[strand]:
                for coord in range(startCoord, endCoord):
                    if coord in self.d[strand][reference]:
                        self.d[strand][reference][coord].append(dataObject)
                    else:
                        self.d[strand][reference][coord] = []
                        self.d[strand][reference][coord].append(dataObject)
            else:
                self.d[strand][reference] = {}
                for coord in range(startCoord, endCoord):
                    self.d[strand][reference][coord] = []
                    self.d[strand][reference][coord].append(dataObject)
        elif self.dataType == list:
            #check that the site has been written to the chromDict in self.d storage object
            if reference in self.d[strand]:
                for coord in range(startCoord, endCoord):
                    if coord in self.d[strand][reference]:
                        value_l = []
                        assert len(self.d[strand][reference][coord]) == len(dataObject)
                        for t in zip(self.d[strand][reference][coord], dataObject):
                            x, y = t
                            if type(y) == float or type(y) == int:
                                assert type(x) == type(y)
                                value_l.append( x + y )
                            elif type(y) == str:
                                assert type(x) == list
                                value_l.append(x.append(y))
                            else:
                                print('bioPyL does not currently support value type:')
                                print(type(y))
                                raise TypeError
                    else:
                        self.d[strand][reference][coord] = []
                        # must use nested list for string types
                        for i in dataObject:
                            if type(i) == int or type(i) == float:
                                self.d[strand][reference][coord].append(i)
                            elif type(i) == str:
                                self.d[strand][reference][coord].append( [i,] )
            else:
                self.d[strand][reference] = {}
                for coord in range(startCoord, endCoord):
                    self.d[strand][reference][coord] = []
                    # must use nested list for string types
                    for i in dataObject:
                        if type(i) == int or type(i) == float:
                            self.d[strand][reference][coord].append(i)
                        elif type(i) == str:
                            self.d[strand][reference][coord].append( [i,] )
        else:
            print('unsupported dataObject dataType!')
            raise IOError
        
        self.write_count += 1
        if self.write_count % self.write_chunk == 0:
            self.d.sync()
            print('{0} sites written'.format(str(self.write_count)))
        
        return


class FastqSQLite():
    """
    SQLite adapters for fastq data.
    """
    
    def __init__(self, fp):
        
        import sqlite3
        
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
    
    def input_fastq(self, fp, col = ('id') ):
        """
        Create a SQLite database from the input fastq file.
        Currently, only the id column is supported.
        """
        #TODO: overwrite existing file, warn user
        #TODO: if useful, implement col = ('id','qual','seq')) 
        #      by putting qual and seq in shelf values.
        #TODO: call get_connection
        
        print('...building fastq database from the file at:')
        print(fp)
        print('...building {0} column(s)...'.format( col ))
        
        from bioPyL.flatfile_parsing import FastqReader
        fq_gen = FastqReader(fp)
        
        if col == ('id'):
            
            def row_gen():
                cnt = 0
                for d in fq_gen:
                    t = (d['TOPID'][1:], ) #strips leading '@'
                    cnt += 1
                    if cnt % 1000000 == 0:
                        print('{0} reads written'.format(str(cnt)))
                    yield t
            
            self.c.execute('''CREATE TABLE reads (
                              read_id TEXT, 
                              PRIMARY KEY(read_id))''')
            
            self.c.executemany('''INSERT INTO reads VALUES 
                                   (?)''', row_gen())
            
            self.conn.commit()
            #TODO: test speed of indexes after the database is full instead of
            # declaring the primary key during db creation
        
        else:
            raise IOError('bioPyL does not yet support that column structure')
        
        return
    
    def id_exists(self, read_id):
        
        # check if the read_id is present
        t = (read_id, )
        self.c.execute('''SELECT EXISTS (
                            SELECT 1 FROM reads WHERE read_id=(?) LIMIT 1
                            )''', (read_id, ))
        
        return bool(self.c.fetchone()[0])

#TODO: try a "build index" routine and see if it is fasteer than the shelf
#...you may be able to standardize on sqlite
class HitsClipSQLite():
    """
    SQLite adapters for fastq data.
    """
    
    def __init__(self, fp):
        
        import sqlite3
        
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
    
    def input_positional_shelf(self, fp):
        """
        Create a SQLite database from the input positional shelf file.
        """
        
        #TODO: overwrite existing file, warn user
        #TODO: check that the input data is integer type
        # only integer type needs to be supported
        #TODO: try running index operation after building the
        #db as opposed to declaring primary key
        
        self.c.execute('''CREATE TABLE clip_signature
                             (
                                ref TEXT, 
                                coord INT, 
                                strand TEXT, 
                                raw_coverage INT, 
                                single_nt_del INT, 
                                cleavage INT, 
                                PRIMARY KEY(ref, coord, strand)
                             )''')
        
        with PositionalShelf(fp) as p_db:
            
            def row_gen():
                print('positive_strand')
                cnt = 0
                strand = '+'
                for ref in p_db.d[strand].keys():
                    for coord in p_db.d[strand][ref].keys():
                        raw_coverage = p_db.d[strand][ref][coord]
                        t = (ref, coord, strand, raw_coverage, 0, 0)
                        cnt += 1
                        if cnt % 1000000 == 0:
                            print('{0} reads written'.format(str(cnt)))
                        yield t
            
            self.c.executemany('''INSERT INTO clip_signature VALUES
                                   (?,?,?,?,?,?)''', row_gen())
            
            def row_gen():
                print('negative_strand')
                cnt = 0
                strand = '-'
                for ref in p_db.d[strand].keys():
                    for coord in p_db.d[strand][ref].keys():
                        raw_coverage = p_db.d[strand][ref][coord]
                        t = (ref, coord, strand, raw_coverage, 0, 0)
                        cnt += 1
                        if cnt % 1000000 == 0:
                            print('{0} reads written'.format(str(cnt)))
                        yield t
            
            self.c.executemany('''INSERT INTO clip_signature VALUES
                                   (?,?,?,?,?,?)''', row_gen())
        
        self.conn.commit()
        
        return
    
    def build_index_from_sam(self, fp):
        """
        Create a SQLite database index from the sam file.
        """
        #SLOW DO NOT USE
        #KEPT FOR POSTERITY TO ACT AS TEMPLATE FOR TRANSCRIPTOME 
        #ALIGNMENT DATA (MANY REFS)
        
        #TODO: overwrite existing file, warn user
        #TODO: check that the input data is integer type
        # only integer type needs to be supported
        #TODO: try running index operation after building the
        #db as opposed to declaring primary key
        
        self.c.execute('''CREATE TABLE clip_signature
                             (
                                ref TEXT, 
                                coord INT, 
                                strand TEXT, 
                                raw_coverage INT, 
                                single_nt_del INT, 
                                cleavage INT
                             )''')
        
        from bioPyL.flatfile_parsing import SamReader
        with SamReader(fp) as sam_gen:
            #TODO: parameterize uniq_only in method call
            def row_gen(uniq_only = True):
                cnt = 0
                for d in sam_gen:
                    if sam_gen.isMapped():
                        #n_of_mapped_reads += 1
                        #identify uniquely aligned reads by accessing 
                        #the sam optdict 
                        if sam_gen.unique_bool():
                            #n_of_unique_aligns += 1
                            pass
                        elif uniq_only:
                            continue
                        else:
                            pass
                        
                        reference, startCoord, endCoord, strand = sam_gen.calc_ucsc_interval()
                        for coord in range(startCoord, endCoord):
                            t = (reference, coord, strand, 0, 0, 0)
                            cnt += 1
                            if cnt % 1000000 == 0:
                                print('{0} reads written'.format(str(cnt)))
                            yield t
            
            self.c.executemany('''INSERT INTO clip_signature VALUES
                                   (?,?,?,?,?,?)''', row_gen())
        #TODO: build index statement...declaring primary key above incurs 
        # need for insert or replace statement. it should be faster to just 
        # index after data insert.
        #PRIMARY KEY(ref, coord, strand)
        
        self.conn.commit()
        
    
    def update_oneD_from_sam(self, fp):
        """
        Update the single_nt_del column using the input sam file.
        """
        
        #TODO: include the call to get_connection from here
        from bioPyL.flatfile_parsing import SamReader
        
        sam_gen = SamReader(fp)
        
        def row_gen():
            cnt = 0
            for d in sam_gen:
                #TODO: make uniq a definable parameter
                if sam_gen.unique_bool():
                    l = sam_gen.get1DsiteList()
                    for siteTup in l:
                        reference, startCoord, endCoord, strand = siteTup
                        t = (reference, startCoord, strand)
                        cnt += 1
                        if cnt % 1000000 == 0:
                            print('{0} reads written'.format(str(cnt)))
                        yield t
        
        #TODO: make this into an update statement
        self.c.executemany('''UPDATE clip_signature SET single_nt_del = single_nt_del + 0 WHERE ref = ? AND coord = ? AND strand = ?''', row_gen())
        
        
        self.conn.commit()
        
        return


###########NOTES
#def compileKnownGeneSQLite():
#    """
#    a function that will compile a bed12 file containing a gene annotation
#    to a sqlite table. note: this function is more generic than indicated by the 
#    function name...beacuse it is not just the knowngene set that can be compiled by this
#    function (any bed12 will work).
#    """
#    
#    #TODO: include all of the bed12 fields ing the compiled database.
#    
#    import sqlite3
#    #http://docs.python.org/3.2/library/sqlite3.html
#    # also see Jay Kreibich's book "Using SQLite" published by O'Reilly
#    
#    from bioPyL.bioFileRW.bedO import Bed12Reader
#    
#    conn = sqlite3.connect(os.path.join(findDataDirectory(), 'knownGene_sqlite.dat'))
#    c = conn.cursor()
#    
#    # Create table
#    c.execute('''CREATE TABLE knownGene
#                 (
#                    ucsc_acc TEXT, 
#                    hg19_chrom TEXT, 
#                    hg19_strand TEXT, 
#                    hg19_start INT, 
#                    hg19_end INT
#                 )''')
#    
#    knownGeneBedFP = findKnownGeneBed12()
#    knowngeneBedgen = Bed12Reader(knownGeneBedFP)
#    
#    for bedEntry in knowngeneBedgen:
#        
#        ucsc_acc = bedEntry['name']
#        hg19_chrom = bedEntry['chrom']
#        hg19_strand = bedEntry['strand']
#        hg19_start = bedEntry['chromStart']
#        hg19_end = bedEntry['chromEnd']
#        
#        rowT  = (ucsc_acc, hg19_chrom, hg19_strand, hg19_start, hg19_end)
#        c.execute('''INSERT INTO knownGene VALUES (?,?,?,?,?)''', rowT)
#    
#    # Save (commit) the changes
#    conn.commit()
#    
#    return


#samFP, uniqOnly = True):
#    print('WARNING: this function inserts the library_stats key into the')
#    print('output dict, therefore you must explicitly declare the')
#    print('strand keys instead of iterating over all top level keys for 
#    print('functions that effect the whole chromDict\n')
#    
#    from bioPyL.bioFileRW.samO import SamReader
#    bwaSamGen = SamReader(samFP)
#    countUniq = 0
#    countAllReads = 0 #note: bwa outputs one line per read (even multialigning reads only get one line by default)
#    countMappedRds = 0
#    chromDict = strandedChromDict()
#    dataObject = 1 # the data object component of the siteTup will always be 1 because we are merely calculating raw coverage
#    
#    for alignD in bwaSamGen:
#        countAllReads += 1
#        
#        if bwaSamGen.isMapped():
#            
#            countMappedRds += 1
#            
#            #identify uniquely aligned reads by accessing the optional sam fields (stored in optdict)
#            uniqBool = bwaSamGen.isUniq()
#            if uniqBool: countUniq += 1
#            
#            if not uniqBool and uniqOnly:
#                pass
#            else:
#                reference, startCoord, endCoord, strand = bwaSamGen.calcUCSCcoordTup()
#                chromDict.addSiteTup( (reference, startCoord, endCoord, strand, dataObject) )
#                #print(reference, startCoord, endCoord, strand, dataObject)#debugging
#    
#    # write the library stats to key
#    chromDict.d['library_stats'] = { 'totalReads' : countAllReads, #TODO: add totalMappedReads field like in bowtie routine
#                                     'totalUniqReads' : countUniq,
#                                     'totalMappedReads' : countMappedRds }
#    
#    return chromDict



#def getGeneNamesForBed6(samFP):
#    
#    """
#    Input a bed6 filepath and a coordinate lookup is performed for each line.
#    The intersecting gene symbols are appended to the name field and a new bed6
#    file is written with .wGeneSymbols extension
#    
#    """
#    print('WARNING: BED FILE MUST BE HEADERLESS!')
#    #TODO: add a header argument that will allow bed files with headers
#    
#    import sqlite3
#    #http://docs.python.org/3.2/library/sqlite3.html
#    # also see Jay Kreibich's book "Using SQLite" published by O'Reilly
#    
#    from bioPyL.data.dataWrapper import findKnownGeneSQLite
#    from bioPyL.data.dataWrapper import ucscAccGeneIDlookupD
#    
#    #####connect to knownGene table db
#    conn = sqlite3.connect(findKnownGeneSQLite())
#    c = conn.cursor()
#    
#    
#    bedGen = Bed6Reader(bedFP)
#    LofLofTup = []
#    outFP = bedFP + '.geneSymbols'
#    
#    for bedEntry in bedGen:
#        t = (bedEntry['chrom'], bedEntry['strand'], bedEntry['chromStart'], bedEntry['chromEnd'])
#        c.execute('''SELECT ucsc_acc FROM knownGene
#                     WHERE hg19_chrom = ?
#                     AND hg19_strand = ?
#                     AND (? BETWEEN hg19_start AND hg19_end OR ? BETWEEN hg19_start AND hg19_end)''', t)
#        LofLofTup.append(c.fetchall())
#    
#    ucscAcclookupD = ucscAccGeneIDlookupD()
#    
#    LofLofSymbols = []
#    for LofTup in LofLofTup:
#        geneSymbolL = []
#        
#        if LofTup == []:
#            geneSymbol = 'UNKNOWN'
#            geneSymbolL.append(geneSymbol)
#        
#        #note: this for loop will not execute in the case of empty list
#        for tup in LofTup:
#            ucscAcc, = tup
#            geneSymbol, geneDsc = ucscAcclookupD.get( ucscAcc, ('UNKNOWN', 'UNKNOWN') )
#            
#            #replace spaces with dashes so that it doesn't screw
#            #the downstream parser
#            #TODO: make this a stronger check for bad characters with regular expressions
#            #http://stackoverflow.com/questions/1007481/how-do-i-replace-whitespaces-with-underscore-and-vice-versa
#            if ' ' in geneSymbol:
#                geneSymbol = geneSymbol.replace(" ","_")
#            
#            geneSymbolL.append(geneSymbol)
#        
#        LofLofSymbols.append(list(set(geneSymbolL)))
#    
#    print(len(LofLofTup)) #debugging
#    print(len(LofLofSymbols)) #debugging
#    
#    #now writeout
#    with open(outFP, 'w') as outFH:
#        
#        bedGen = Bed6Reader(bedFP)
#        for bedEntry, LofSymbols in zip(bedGen,LofLofSymbols):
#            t = ( bedEntry['chrom'], 
#                  str(bedEntry['chromStart']), 
#                  str(bedEntry['chromEnd']), 
#                  bedEntry['name'] + '_' + '_'.join([ geneSymbol for geneSymbol in LofSymbols]), # append gene symbols
#                  str(bedEntry['score']), 
#                  bedEntry['strand'] )
#            
#            outFH.write( '\t'.join(t) + '\n')
#    
#    # We can also close the cursor if we are done with it
#    c.close()
#    
#    return


