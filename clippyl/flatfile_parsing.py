# -*- coding: utf-8 -*-
#!usr/bin/python3

import re
#from bioPyL.bioUtils.compareCoords import getOverlap
#from bioPyL.bioFileRW.bedO import bed12_reader

#TODO: FIX FASTQREADER __init__ so that the program
# exits gracefully if a bad fp is passed
class FastqReader():
    """
    This iterable object class will read a fastq file entrywise
    and yield a dictionary containing read entry info.
    """
    
    # USE FOR 'SINGLE LINE' FASTQ ONLY
    # the FastqReader expects 4 lines per read
    # hint: if you do not have a 'single line' fastq file then 
    # you may need to preprocess with fastx-toolkit's fasta_formatter
    
    #TODO: write a method to detect quality encoding
    #TODO: write a method to decode quality encoding
    
    def __init__(self, filepath):
        self.header = []
        try:
            self.fh = open(filepath)
            # read each line of read entry
            self.topid = self.fh.readline().rstrip() # ID entry
            self.read = self.fh.readline().rstrip() # read
            self.botid = self.fh.readline().rstrip() # ID entry
            self.qual = self.fh.readline().rstrip() # quality string
        except IOError as ioerr:
            "{0} has no data or does not exist".format(filepath)
    
    def __enter__(self):
        return self
    
    def __exit__(self, *exc_info):
        self.fh.close()
    
    def __iter__(self):
        return self
    
    def __next__(self):
        if self.read:
            self.d = {'TOPID' : self.topid, # ID entry
                      'READ' : self.read, # read
                      'BOTID' : self.botid, # ID entry
                      'QUAL' : self.qual # quality string
                      }
            self.topid = self.fh.readline().rstrip() # ID entry
            self.read = self.fh.readline().rstrip() # read
            self.botid = self.fh.readline().rstrip() # ID entry
            self.qual = self.fh.readline().rstrip() # quality string
            return self.d
        else:
            raise StopIteration
    
    def format_validator(self):
        '''
        Returns false if the fastq header flags are not as expected.
        '''
        if self.d['TOPID'][0] != '@' or self.d['BOTID'][0] != '+':
            return False
        else:
            return True
    
    def detect_pair_info(self):
        '''
        The read ids sometimes end with /1 or /2 to denote paired-end or
        mate-pair reads. it can be usefule to remove it. 
        ex: bwa seems to remove this tag from the ids in the 
        output sam file that it generates so I need to remove it to compare
        ids with the fastq file it was derived from.
        '''
        #http://en.wikipedia.org/wiki/FASTQ_format
        
        if self.d['TOPID'][-2] == '/':
            return True
        else:
            return False
    
    def remove_pair_info(self):
        '''
        The read ids sometimes end with /1 or /2 to denote paired-end or
        mate-pair reads. You can look for it with detect_fq_pair_info. 
        usage example: bwa seems to remove this tag from the ids in the 
        output sam file that it generates so I need to remove it to compare
        ids with the fastq file it was derived from.
        '''
        #http://en.wikipedia.org/wiki/FASTQ_format
        self.d['TOPID'] = self.d['TOPID'][:-2]
        # note that BOTID is an optional field. If it is not present then
        # it can result in index error here when sliced.
        try:
            self.d['BOTID'] = self.d['BOTID'][:-2]
        except IndexError:
            pass
        
        return
    
    def validate_readid(self):
        '''Given a read id with a valid format, the format type will 
           be returned. An error will be produced if it is not a valid 
           read_id type'''
        
        # The parsing logic is based on the SolexaQA code.
        # corresponds to ./LB_parse_qname_v2.3.pl
        # http://solexaqa.sourceforge.net/
        # NOTE: LB user group post: 
        # https://groups.google.com/d/topic/solexaqa-users/Wj_fZbBqK2o/discussion
        # FAQ: http://solexaqa.sourceforge.net/questions.htm#headers
        # Translated to python:
        # Perl regex help: http://www.cs.tut.fi/~jkorpela/perl/regexp.html
        # https://docs.python.org/3.2/howto/regex.html
        # All possible fields are indicated at the wikipedia page.
        #TODO: get definition document from Illumina
        #http://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers
        # set an error message to use when malformed entries are encountered
        error_message = '''CLIP-PyL does not recognize the read id format.\n
                           The read id that CLIP-PyL attempted to parse is:\n
                           {0}\n'''
        
        tile_number = 0
        space_split = self.d['TOPID'].split(' ')
        if re.match(r'^\@[ES]R[RA][\d]+.[\d]+\s*', self.d['TOPID']) and len(space_split) > 1:
            l = space_split[1].split(':')
            if len(l) < 3:
                raise IOError(error_message.format(self.d['TOPID']))
            if len(l) < 8:
                tile_number = l[-3]
                format_type = 'typeA'
            else:
                tile_number = l[-4]
                format_type = 'typeB'
        elif self.d['TOPID'][0] == '@':
            l = space_split[0].split(':')
            if len(l) < 3:
                raise IOError(error_message.format(self.d['TOPID']))
            if len(l) < 8:
                tile_number = l[-3]
                format_type = 'typeC'
            else:
                tile_number = l[-4]
                format_type = 'typeD'
        else:
            raise IOError(error_message.format(self.d['TOPID']))
        
        if tile_number == 0:
            raise IOError(error_message.format(self.d['TOPID']))
        else:
            print('Read ID format detected: {0}'.format(format_type))
        
        return format_type

def validate_fastq_format(fp):
    with FastqReader(fp) as g:
        print('Validating format for file:')
        print(fp)
        g.__next__()
        if g.format_validator():
            print('fastq format recognized')
            g.validate_readid()
            return True
        else:
            return False

def detect_fq_pair_info(fp):
    '''
    The read ids sometimes end with /1 or /2 to denote paired-end or
    mate-pair reads. It can be useful to remove it. 
    ex: bwa seems to remove this tag from the ids in the 
    output sam file that it generates so I need to remove it to compare
    ids with the fastq file it was derived from.
    '''
    #http://en.wikipedia.org/wiki/FASTQ_format
    with FastqReader(fp) as g:
        print('Checking for matepair tag:')
        print(fp)
        g.__next__()
        if g.detect_pair_info():
            print('mate pair tag detected')
            return True
        else:
            return False

#TODO:Rename this class to bwa-sam? because all method work on BWA-derived Sam files.
class SamReader():
    """
    This function will read a sam file
    by creating generator object in memory.
    The sam file standard is here:
    http://samtools.sourceforge.net/samtools.shtml
    http://samtools.sourceforge.net/SAM1.pdf
    Sam file standard adhered to by bwa is here:
    http://bio-bwa.sourceforge.net/bwa.shtml
    
    WARNING: OPTIONAL FIELDS REPORTED DIFFER ACROSS ALIGNER THAT GENERATED THE FILE
    WARNING: STYLE OF MULTI-READ REPORTING DIFFERS ACROSS ALIGNER AS WELL
    
    provide path to a bed12 file as overlapCheck arg and you can call
    the method queryLookUpDict to see if there is overlap between the current
    self.d['POS'] and any entry in the bed12 file
    """
    #TODO: Parse @PG header line to define optional field behavior
    
    def __init__(self, filepath, overlapCheck = None):
        self.header = []
        try:
            self.fh = open(filepath)
            self.line = self.fh.readline()
            while self.line.find('@') == 0:
                #TODO: parse headr to nested dicts with @ tag and @SEQ lines to 
                # a dict object with chrom names as keys and lengths (int), 
                self.header.append(self.line)
                self.line = self.fh.readline()
            self.alignment = self.line
        except IOError as ioerr:
            print('{} not found'.format(filepath))
        if overlapCheck:
            print('Overlap check mode activated')
            print('building lookup structure for file:\n {0}'.format(overlapCheck))
            self.lookUpDict = self.buildLookUpFromBed(overlapCheck)
            print('you can now use the queryLookUpDict method')
        ###
        #CIGAR string parsing regex
        ###
        # I hard coded max read length into the regex pattern.
        # for read lengths longer than 39 bp you must alter the regex
        # also, the bitwise flag decoding is more complicated for paired end reads
        # (see strand determination code below; i['FLAG']
        #TODO: Use [0-9]+ instead of ([0-9]|[0-9][0-9]) and you can remove this warning
        print('WARNING: USE ONLY FOR SINGLE-END <=99 BP SINGLE END READS')
        # see special characters '[',']' and '|' http://docs.python.org/library/re.html#regular-expression-syntax
        # http://docs.python.org/dev/howto/regex.html
        # http://docs.python.org/howto/regex.html#grouping
        self.CIGARregexP = re.compile(r'([0-9]|[0-9][0-9])([MIDNSHP=X])')
        #minimal example
        #CIGAR = '14M1D22M'
        #m = CIGARregexP.findall(CIGAR)
        #m
        #for count, operation in m:
        #    foo
        
        #TODO: create an MD:Z mismatch string parsing Regex
        #[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*
        #The MD field aims to achieve SNP/indel calling without looking at the reference. For example, a string `10A5^AC6'
        #means from the leftmost reference base in the alignment, there are 10 matches followed by an A on the reference which
        #is different from the aligned read base; the next 5 reference bases are matches followed by a 2bp deletion from the
        #reference; the deleted sequence is AC; the last 6 bases are matches. The MD field ought to match the CIGAR string.
        self.MDZregexP = re.compile(r'[0-9]+|[A-Z]+|\^[A-Z]+')
        
    
    def __enter__(self):
        return self
    
    def __exit__(self, *exc_info):
        self.fh.close()
    
    def __iter__(self):
        return self #TODO: can you just put a yield statement in the __next__ method?
    
    def __next__(self):
        if self.alignment:
            self.raw = self.alignment
            self.l = self.alignment.rstrip().split()
            self.d = {'QNAME' : self.l[0], # Name of read
                      'FLAG' : int(self.l[1]), # 
                      'RNAME' : self.l[2], # reference sequence name (e.g. chr11)
                      'POS' : int(self.l[3]), # 
                      'MAPQ' : int(self.l[4]), # 
                      'CIGAR' : self.l[5], # 
                      'RNEXT' : self.l[6], # 
                      'PNEXT' : int(self.l[7]),
                      'TLEN' : int(self.l[8]),
                      'SEQ' : self.l[9], #NOTE: The seq that results from all CIGAR string operations is shown. Also, if the read maps to the reverse-strand then the reverse complement of the read is shown
                      'QUAL' : self.l[10],
                      'OPTDICT' : {}, #instantiate dict for optional fields; note that an empty dict will return
                      'OPTLIST' : '' #instatiate list for optionalfields, facillitates easy writeout
                      }
            # Build a nested dict from optional fields, if they are present
            if len(self.l) > 11:
                self.optDict = {}
                self.d['OPTLIST'] = self.l[11:]
                for optField in self.l[11:]:
                    optList = optField.split(':')
                    self.d['OPTDICT'][optList[0]] = (optList[1], optList[2]) # creates a nested dictionary of tuples, the dict is hashed by opt field name
            # Add strand info
            self.d['STRAND'] = self.gleanStrand()
            # load the next line for next iteration
        else:
            print('done')
            raise StopIteration
        try:
            self.alignment = self.fh.readline()
        except ValueError:
            self.alignment = None
            return self.d
        return self.d
    
    def headerOut(self):
        #format header line list into a string for easy text file writout
        headerStr = ''
        for headerLine in self.header:
            headerStr += headerLine
        return headerStr
    
    def samLineOut(self):
        #return current sam alignment string for easy text file writout
        return self.raw
    
    def dictHead(self):
        dictHead = {} # we will populate this dict with dicts to map parsed header info
        for headerLine in self.header:
            l = headerLine.split()
            recType = l[0][1:] # parse record type code
            if recType not in dictHead: dictHead[recType] = {}
            tagValLisOfTups = [tuple(tagValpair.split(':')) for tagValpair in l[1:]] # parse tag:value pairs
            
            if recType == 'SQ': # record type is reference sequence dictionary
                for tag, value in tagValLisOfTups:
                    if tag == 'SN':
                        refSeqName = value
                    elif tag == 'LN':
                        refSeqLen = value
                    else:
                        pass # discard other fields
                #TODO: use orderedDict for this and preserve order which defines alignment sort order
                dictHead['SQ'][refSeqName] = refSeqLen
            else:
                pass #TODO:add elif statements to catch other recTypes
        return dictHead
    
    def gleanStrand(self):
        #'WARNING: COMPATIBLE WITH SINGLE END READS ONLY!'
        # make strand-aware bedgraph files by decoding bitwise flag and appending strand info to the site entry
        # top strand will be written to fhOutBg1Dplus and bottom strand to fhOutBg1Dminus.
        # sam bit flag http://seqanswers.com/forums/showthread.php?t=17314
        # TODO: make a generalized function for converting decimal to 
        # binary and generate a reference dictonary to lookup code definitions
        # decoding bitwise flags: http://seqanswers.com/forums/showthread.php?t=17314
        # and http://seqanswers.com/forums/archive/index.php/t-18263.html
        # and http://biostar.stackexchange.com/questions/13079/flag-in-sam-format
        if self.d['FLAG'] == 0 or self.d['FLAG'] == 256: strand = '+'
        elif self.d['FLAG'] == 16 or self.d['FLAG'] == 272: strand = '-'
        elif self.d['FLAG'] == 4 or self.d['FLAG'] == 20: strand = 'NA' #unmapped read
        else: raise NameError('Did not understand flagbit {0}'.format(str(self.d['FLAG'])))
        return strand
    
    def isMapped(self):
        #'WARNING: COMPATIBLE WITH SINGLE END READS ONLY!'
        # make strand-aware bedgraph files by decoding bitwise flag and appending strand info to the site entry
        # top strand will be written to fhOutBg1Dplus and bottom strand to fhOutBg1Dminus.
        # sam bit flag http://seqanswers.com/forums/showthread.php?t=17314
        # TODO: make a generalized function for converting decimal to 
        # binary and generate a reference dictonary to lookup code definitions
        # decoding bitwise flags: http://seqanswers.com/forums/showthread.php?t=17314
        # and http://seqanswers.com/forums/archive/index.php/t-18263.html
        # and http://biostar.stackexchange.com/questions/13079/flag-in-sam-format
        
        if self.d['FLAG'] == 4 or self.d['FLAG'] == 20: #unmapped read
            isMappedBool = False
        else: 
            isMappedBool = True
        return isMappedBool
    
    def parseCIGAR(self):
        
        # parse the CIGAR string for the current read, return list of tuples
        # comprised of positionCount, operation
        # WARNING: this returns the native CIGAR one-based coord system
        m = self.CIGARregexP.findall(self.d['CIGAR'])
        return m
    
    def parseMDZ(self):
        
        #parse the MD:Z field containing the mismatch string
        #NOTE: you must check to be sure the read is mapped before calling
        #this method because you will get a key error if there is no
        #optdict, like when the read is unmapped
        sam_data_type, mismatch_str = self.d['OPTDICT']['MD']
        m = self.MDZregexP.findall(mismatch_str)
        return m
    
    def detect_mm_variant(self, native_nt = 'T', variant_nt = 'C', 
                                siteTup_mode = False, 
                                nt_set = 'ATCGN'):
        '''
        you can check a specific nt variant.
        the default mismatch that is detected is T to C.
        returns a list of integers for each mismatch event.
        the integers returned are the zero-based indices of the variant 
        occurence along the string stored in SEQ field.
        a list of zero-based chromosomal coordinates can be returned instead 
        by setting siteTup_mode to True. the chromosomal coordinates will
        be returned as a list of siteTups can be calculated.
        '''
        v = []
        #NOTE: this read must be mapped. otherwise optdict will be empty.
        sam_data_type, mismatch_str = self.d['OPTDICT']['MD']
        if native_nt in mismatch_str:
            m = self.MDZregexP.findall(mismatch_str)
            sam_seq_pos = -1 # this is the zero-based index into the SEQ string
            ref_pos = self.d['POS'] - 2 # this creates the zero-based index into the REF
            for i in m:
                if native_nt in i and '^' != i[0]:
                    for nt in i:
                        sam_seq_pos += 1
                        ref_pos += 1
                        if nt == native_nt:
                            if self.d['SEQ'][sam_seq_pos] == variant_nt:
                                #print('found one', self.d['QNAME']) #debugging
                                if siteTup_mode:
                                    v.append( (self.d['RNAME'], ref_pos, ref_pos + 1, self.gleanStrand()) )
                                else:
                                    v.append( (sam_seq_pos, ref_pos) )
                            else:
                                continue
                        else:
                            continue
                
                elif '^' == i[0]:
                    ref_pos -= len(i[1:])
                    continue
                
                elif i in nt_set:
                    sam_seq_pos += len(i)
                    ref_pos += len(i)
                    continue
                
                elif re.match('[0-9]+', i):
                    sam_seq_pos += int(i)
                    ref_pos += int(i)
                    continue
                
                else:
                    print('something went wrong.')
                    raise IOError
        
        return v
    
    def calc_ucsc_interval(self, cleavage_coords = False):
        """
        a method to calculate UCSC-style zero based, right open coords from the one-based SAM coord,
        it returns (chrom, startCoord, endCoord, strand)
        
        optionally, return a tuple of two tuples mapping the fiveprime and threeprime terminal nucleotides
        with the cleavage_coords argument
        
        """
        reference = self.d['RNAME']
        # glean UCSC-style coord aka (zero-based, right open]
        startCoord = self.d['POS'] - 1 #calculate UCSC-style zero-based start coordinate 
        # must access the CIGAR operation to calclulate fragment length. otherwise we would not account for
        # deletions in the length of the actual fragment
        m = self.parseCIGAR()
        # get sum of CIGAR operations posCounts
        # (excluding insertion operations that do not refernce the read index; they reference the reference 
        # index). This will yield the index (one-based) for the site of the deleted base in the 
        # read string.
        # see get1DsiteList method below for logic to extract coords of nucleotides that underwent specific mapping operations
        # (e.g. deletion, insertion)
        realFragLength = sum([int(posCount) for posCount, operation in m if operation != 'I'])
        endCoord = startCoord + realFragLength #calculate UCSC-style open right coordinate
        strand = self.gleanStrand()
        
        if cleavage_coords: # return a tuple of siteTuples of the terminal nucleotides for the current fragment
            # first, a siteTup of the 5' terminal nucleotide of the fragment
            # second, a siteTup of the 3' terminal nucleotide of the fragment
            if strand == '+':
                siteTup5pr = ( reference, startCoord, startCoord + 1, strand )
                siteTup3pr = ( reference, endCoord - 1, endCoord, strand )
            elif strand == '-':
                siteTup5pr = ( reference, endCoord - 1, endCoord, strand )
                siteTup3pr = ( reference, startCoord, startCoord + 1, strand )
            else:
                print('unrecognized strand...quitting')
                print(strand)
                raise IOError
            return( siteTup5pr, siteTup3pr )
        
        else: #default behavior, return a single siteTup
            return ( reference, startCoord, endCoord, strand )
    
    def get1DsiteList(self):
        # return a list of tuples of single basepair deletion sites
        # gleaned from the current alignment line
        # with chrom, startCoord, endCoord, strand
        # note: coords are converted and zero-based coords are returned
        oneDSitesList = []
        m = self.parseCIGAR()
        # http://stackoverflow.com/questions/522563/accessing-the-index-in-python-for-loops
        for idx, tup in enumerate(m):
            posCount, operation = tup
            if posCount == 1 and operation == 'D':
                # get the sum of posCounts for all CIGAR operations previous to the
                # found deletion (excluding insertion operations that do
                # not refernce the read index; they reference the reference
                # index). This will yield the index (one-based) for the site of the deleted base in the
                # read string.
                #readPos = sum([int(tup[0]) for tup in m[:idx] if tup[1] != 'D'])
                # same as above, but includes insertion operations to yield the index (aka coordinate)
                # of the deleted base in the reference genome (one-based).
                refPosOffset = sum([int(t[0]) for t in m[:idx + 1] if t[1] != 'I'])
                delCoordZeroBase = self.d['POS'] + refPosOffset - 2 # deduce zero-based coordinate for bedgraph writeout
                #note 2 is subtrated because both the start pos and the offset are onebased (one for each)
                # end-coord is lazily claculated by simply adding one, this works for single base deletions but
                # a function that is generalixable to multiple base sites would use the CIGAR code instead
                oneDSitesList.append((self.d['RNAME'], delCoordZeroBase, delCoordZeroBase + 1, self.gleanStrand()))
                # NOTE: The CIGAR operation M character does not encode mismatch
                # you must consult the MD tag for this info. For example:
                # 14M1D22M ATGATAACTAATGACAAAAAAAAAAAAAAAAAAAAA
                # means deletion at 14 ... calc position with 
                # MD:Z:9G4^G22
                # http://seqanswers.com/forums/showthread.php?t=2174
                #TODO: report deleted base, compare most commonly deleted base with other datasets to glean x-link nt bias
                #NOTE: the deleted base can be gleaned by looking at the MD field
                # so you don't have to do a lengthy refernce lookup
                # MD field aims to achieve SNP/indel calling without looking at
                # the reference. For example, a string ‘10A5^AC6’ means from the
                # leftmost reference base in the alignment, there are 10 matches
                # followed by an A on the reference which is different from the
                # aligned read base; the next 5 reference bases are matches 
                # followed by a 2bp deletion from the reference;
                # the deleted sequence is AC;
                # the last 6 bases are matches
        return oneDSitesList
        
    def get1DoffsetList(self):
        # for use in Chaolin's permutation test
        # return a list of integers, each representing one single basepair deletion site
        # gleaned from the current alignment line
        # note: offsets are converted from SAM's native one-based index and zero-based offsets are returned
        oneDoffsetList = []
        m = self.parseCIGAR()
        # http://stackoverflow.com/questions/522563/accessing-the-index-in-python-for-loops
        for idx, tup in enumerate(m):
            posCount, operation = tup
            if posCount == 1 and operation == 'D':
                # get the sum of posCounts for all CIGAR operations previous to the
                # found deletion (excluding insertion operations that do
                # not refernce the read index; they reference the reference
                # index). This will yield the index (one-based) for the site of the deleted base in the
                # read string (aka offset)
                readPos = sum([int(tup[0]) for tup in m[:idx+1] if tup[1] != 'D'])
                # convert to zero based and append to list
                oneDoffsetList.append(readPos-1)
        return oneDoffsetList
    
    def unique_bool(self, samType = 'bwa'): #TODO: move the samType argument to the init routine, parse it from the headers
        """report boolean, give true if uniquely aligned read"""
        #determine unique alignment by accessing the optional sam fields (stored in optdict)
        #(note for samType bwa, the code is reused from samUniqifier function below)
        
        b = 0
        if samType == 'bwa':
            if 'XT' in self.d['OPTDICT']:
                if self.d['OPTDICT']['XT'][1] == 'U':
                    b = 1
        elif samType == 'bowtie':
            print('Bowtie does not write XT optional flag field, the only way to identify uniq reads in')
            print('the bowtie-outputted sam is to sort the file by read ID and check for multiple alihnment entries')
            print('for the same read ID.')
            # see /storage/Ziggy_BigGuy/LB_Bioinformatics/LB_Scripts/seq_data_uniqified_sam_SLBP_CLIP
            # for the old routine that I wrote...I also copied it to ../BioPyL-0.1dev/incorporateme/LB_bowtie_sam_uniqifier_v00
            raise IOError
        else:
            print('Input file format not supported yet, sorry')
            raise IOError
        
        return b
    
    def queryLookUpDict(self):
        #TODO: build this method into the chromDictO data class
        ###############################################################################
        # searchs lookup dict for overlap
        # 
        strand = self.gleanStrand()
        matchBool = 0
        if self.d['RNAME'] in self.lookUpDict[strand]:
            for bedStart, bedEnd in self.lookUpDict[strand][self.d['RNAME']]:
                readStart = self.d['POS'] -1 # convert one-based coord to zero-based
                readEnd = readStart + len(self.d['SEQ']) # right closed end coord
                if getOverlap([readStart, readEnd],[bedStart,bedEnd]) != 0:
                        matchBool = 1
                        break
        return matchBool

def sam_hitsclip_gen(fp, uniq_only = True):
    
    sam_gen = SamReader(fp)
    
    for d in sam_gen:
        
        if sam_gen.isMapped():
            
            #identify uniquely aligned reads by accessing 
            #the sam optdict 
            if sam_gen.unique_bool():
                pass
            elif uniq_only:
                continue
            else:
                pass
            
            #TODO: increment raw_coverage for each nt in the aligned interval
            reference, startCoord, endCoord, strand = sam_gen.calc_ucsc_interval()
            for i in range(startCoord, endCoord):
                #check that the nt exists
                t = (str(reference), i, str(strand), str(reference), int(i), str(strand))
                yield t


#TODO: include method to input a file with the fastq ids only in the
# FastqSQLite class
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
                              read_id TEXT not null, 
                              PRIMARY KEY(read_id))''')
            
            self.c.executemany('''INSERT INTO reads VALUES 
                                   (?)''', row_gen())
            
            self.conn.commit()
            #TODO: test speed of indexes after the database is full instead of
            # declaring the primary key during db creation
        
        else:
            raise IOError('clipPyl does not yet support that column structure')
        
        return
    
    def id_exists(self, read_id):
        
        # check if the read_id is present
        t = (read_id, )
        self.c.execute('''SELECT EXISTS (
                            SELECT 1 FROM reads WHERE read_id=(?) LIMIT 1
                            )''', (read_id, ))
        
        return bool(self.c.fetchone()[0])

class Bed6Reader():
    """
    This class creates a generator object that
    yields a dictionary for each line in a bed file
    (e.g. bedfile representation of UCSC knowngene table).
    The bed file standard is defined by UCSC.
    http://genome.ucsc.edu/FAQ/FAQformat.html#format1
    Each bed line describes a gene model.
    see listserv correspondence for details of annotation coordinate system.
    trackline description is here:
    http://genome.ucsc.edu/goldenPath/help/customTrack.html#TRACK
    """
    
    def __init__(self, inBedFP):
        print("WARNING: REMOVE HEADER FROM BEDFILE BEFORE USING THIS CLASS!")
        print("WARNING: BEDFILE HEADER IS NOT PARSED")
        #TODO: Detect header by looking for word browser or track or # at 
        # beginning of first lines. see:
        # http://genome.ucsc.edu/FAQ/FAQformat.html#format1
        # http://bedtools.readthedocs.org/en/latest/content/overview.html#headers-are-allowed-in-gff-and-bed-files
        
        try:
            self.fh = open(inBedFP, 'r')
            self.bedline = self.fh.readline()
        except IOError as ioerr:
            print("{} not found".format(filepath))
    
    def __enter__(self):
        return self
    
    def __exit__(self, *exc_info):
        self.fh.close()
    
    def __iter__(self):
        return self
    
    def __next__(self):
        if self.bedline:
            self.l = self.bedline.rstrip().split('\t')
            self.d = {'chrom' : self.l[0],
                      'chromStart': int(self.l[1]), #start coord of transcribed fragment for + 
                      'chromEnd' : int(self.l[2]),
                      'name' : str(self.l[3]),
                      'score' : float(self.l[4]), #TODO: may need int() in some instances ?!
                      'strand' : self.l[5]}
        else:
            print('done')
            raise StopIteration
        try:
            self.bedline = self.fh.readline() # load the next line for next iteration
        except ValueError:
            self.bedline = None
            return self.d
        return self.d
    
    #TODO: test this code!
    def calcGraphParams(self, flankint = 0):
        """
        calculate the parameters needed for a graph of the genemodel and surrounding region (specified by flankint).
        """
        # find template code at ./BioPyL/BioPyL-0.1dev/incorporateme/LB_Plotter_Master_v16/LB_code_library/01_Gene_Model_Data_Constructor_v02.py
        
        gmIntervalStart = self.d['chromStart'] - flankint
        gmIntervalEnd = self.d['chromEnd'] + flankint
            
        # deriving important parameters for gene model construction
        # According to the bed12 file format we know that the fields represent the following values:
        # chromosome, chromStart, chromEnd, name, score, strand, thickStart (start codon), thickEnd (stop codon), itemRgb, blockCount (number of exons), blockSizes (exon sizes), blockStarts (exons starts)
        
        priTranscriptLength = self.d['chromEnd'] - self.d['chromStart']
        
        # NOTE: this has some redundant entries 
        # NOTE: I have denoted "absolute for entries that are in absolute chromosome coords.
        dataDict = {'name' : self.d['name'],
                    'chrom' : self.d['chrom'], #absolute
                    'chromStart' : self.d['chromStart'], #absolute
                    'gmIntervalStart' : gmIntervalStart, #absolute
                    'chromEnd' : self.d['chromEnd'], #absolute
                    'gmIntervalEnd' : gmIntervalEnd, #absolute
                    'strand' : self.d['strand'], #absolute
                    }
        
        return dataDict


if __name__ == "__main__":
    pass
