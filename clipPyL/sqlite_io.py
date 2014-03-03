#!/usr/bin/python3
# -*- coding: utf-8 -*-

import sqlite3

#cims_db.py
#cleavage_db.py
#raw_coverage_db.py

#NOTE: this function ONLY includes uniquely aligning reads
# this function computes raw coverage
#TODO: add SamGen.isMapped() logic and totalMappedReads dict key as 
#shown in bowtieSam2rawCoverCD
#TODO: detect sam file type with SamReader method, invoke proper parse scheme


class Fastq2sqlite():
    """
    This class defines methods for interacting with bioPyL's positional 
    database tables. The positional database columns are: 
    reference, coordinate, strand, datum1...datumN
    Thus, these tables are designed to hold data for each position across
    some reference interval (e.g. chromosome, gene model, protein ...).
    """
    
    def __init__(self):
        pass
    
    def get_connection(self, fp):
        
        self.conn = None
        self.conn = sqlite3.connect(fp)
        self.c = self.conn.cursor()
        self.c.execute('SELECT SQLITE_VERSION()')
        self.data = self.c.fetchone()
        print('Loaded SQLite database file at:')
        print(fp)
        print('SQLite version: {0}'.format(str((self.data))))
        #TODO: print tables with # cols and rows
        return
    
    def sam_2_hitsclip_db(self, sam_fp, 
                            cleav_fp, 
                            cleav_ft = 'fq', 
                            uniq_only = True,
                            chunk_size = 100000):
        """
        Create a sqlite HITS-CLIP database.
        Raw coverage and single nucleotide deletion coverage will be gleaned 
        from the input sam file. Cleavage information will be gleaned from the 
        list of adapter-clipped reads.
        NOTE: coordinate system is zero-based
        """
        #TODO: check for file and warn of overwrite (use get_connection)
        
        #instantiate counters for alignment stats table
        n_of_reads = 0 #total number of reads is gleaned from sam file
        n_of_cleaved_reads = 0 #gleaned from cleav_fp
        n_of_mapped_reads = 0 #gleaned from sam flagbits
        n_of_unique_aligns = 0 #gleaned from sam XT opt field
        
        print('Building a CLIP database from the input file at:')
        print(sam_fp)
        from bioPyL.flatfile_parsing import sam_hitsclip_gen
        sam_gen = sam_hitsclip_gen(sam_fp)
        
#        if cleav_ft == 'fq':
#            print('Cleavage info will be gleaned from the .fastq file')
#            print('containing the ADAPTER-CLIPPED READS ONLY at:')
#            print(cleav_fp)
#            from bioPyL.flatfile_parsing import FastqReader
#            fq_gen = FastqReader(cleav_fp)
#            # Create table
#            self.c.execute('''CREATE TABLE cleaved_reads (
#                              read_id TEXT, 
#                              PRIMARY KEY(read_id))''')
#            
#            for d in fq_gen:
#                n_of_cleaved_reads += 1
#                t = (d['TOPID'][1:], )
#                self.c.execute('''INSERT INTO cleaved_reads VALUES (?)''', t)
#            
#            #self.c.execute('''CREATE INDEX idx_cleaved_reads_read_id ON cleaved_reads (read_id)''')
#            # Save (commit) the changes
#            self.conn.commit()
#            
#        
##        elif cleav_ft == 'id':
##            print('Cleavage info will be gleaned from the text file')
##            print('containing the IDs of the ADAPTER-CLIPPED READS ONLY at:')
##            print(cleav_fp)
##            #TODO: implement me
#        
#        else:
#            raise IOError('bioPyL does not support {0} files'.format(cleav_ft))
        
        # Create table
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
        
        self.c.executemany("""INSERT OR REPLACE INTO clip_signature ( ref, coord, strand, raw_coverage, single_nt_del, cleavage) VALUES ( ?, ?, ?, COALESCE((SELECT raw_coverage FROM clip_signature WHERE ref = ? AND coord = ? AND strand = ?), 0), 0, 0)""", sam_gen)
        
        #TODO: move these to __exit__ method
        self.conn.commit()
        self.conn.close()
                
#                #check if contains 1D
#                #increment 1D nucleotides in 1D column
#                
#                
#                #check if it is a cleaved_read
#                t = (d['QNAME'], )
#                self.c.execute('''SELECT EXISTS(SELECT 1 FROM cleaved_reads WHERE read_id=(?) LIMIT 1)''', t)
#                if bool(self.c.fetchone()[0]):
#                    #increment cleav column +1 for terminal nucleotides
#                    siteTup5pr, siteTup3pr = sam_gen.calc_ucsc_interval(cleavage_coords = True)
#                    reference, startCoord, endCoord, strand = siteTup5pr
#                    #increment 5 prime cleavage site
#                    
#                    reference, startCoord, endCoord, strand = siteTup3pr
#                    #increment 3 prime cleavage site
#                    
#                    
        
        
        
        
        #TODO: insert explicit close connection statement!
        
        return



#    def increment_nt_data(self, siteTup):
#        """
#        a method for adding data from a tuple
#        by default the addition operation is used for dataType int or float
#        by default a dataType of str causes the string objects to be aggregated into lists at each site
#        """
#        
#        reference, startCoord, endCoord, strand, dataObject = siteTup
#        #check that the dataObject is expected type
#        if type(dataObject) != self.dataType:
#            print('incorrect data type given in dataObject component of the siteTup')
#            raise IOError
#        
#        #perform the correct data operation for the type given!
#        if self.dataType == int or self.dataType == float:
#            #check that the site has been written to the chromDict in self.d storage object
#            if reference in self.d[strand]:
#                for coord in range(startCoord, endCoord):
#                    if coord in self.d[strand][reference]:
#                        #TODO: use a decorator function here instead that will invoke the correct data 
#                        # object formatting routine for the dataType declared (eg this may call for an append if there are labels stored as a list)
#                        #this would allow you to move the if self.dataType == int: condition elsewhere ?!
#                        self.d[strand][reference][coord] += dataObject
#                    else:
#                        self.d[strand][reference][coord] = dataObject
#            else:
#                self.d[strand][reference] = {}
#                for coord in range(startCoord, endCoord):
#                    self.d[strand][reference][coord] = dataObject
#        elif self.dataType == str:
#            #check that the site has been written to the chromDict in self.d storage object
#            if reference in self.d[strand]:
#                for coord in range(startCoord, endCoord):
#                    if coord in self.d[strand][reference]:
#                        #TODO: use a decorator function here instead that will invoke the correct data 
#                        # object formatting routine for the dataType declared (eg this may call for an append if there are labels stored as a list)
#                        #this would allow you to move the if self.dataType == int: condition elsewhere ?!
#                        self.d[strand][reference][coord].append(dataObject)
#                    else:
#                        self.d[strand][reference][coord] = []
#                        self.d[strand][reference][coord].append(dataObject)
#            else:
#                self.d[strand][reference] = {}
#                for coord in range(startCoord, endCoord):
#                    self.d[strand][reference][coord] = []
#                    self.d[strand][reference][coord].append(dataObject)
#        else:
#            print('unsupported dataObject dataType!')
#            raise IOError
#        
#        return



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


