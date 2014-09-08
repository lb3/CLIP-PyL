import pysam
from clippyl.ome_dict_io import SiteTuple, OmeDict
from clippyl.pysam_cbs import (fetch_hitsclip_dict, 
                               fetch_iclip_dict, 
                               fetch_parclip_dict)

class build_hitsclip_vectors():
    '''
    Generate basewise counts for 1D, cleavage and raw coverage from the 
    set of alignments from the reference interval defined in the site_tup.
    This function parses the list of dictionaries generated by the pysam
    indexed bam file interface function at ??.
    Note: if the cleav_db_fp is set to None then it is assumed that all
    the reads in the sam file are adapter-clipped
    '''
    
    def __init__(self):
        
        self.stat_dict = {}
        self.stat_dict['n_of_reads'] = 0
        self.stat_dict['n_of_strand_matched_reads'] = 0
        self.stat_dict['n_of_unique_aligns'] = 0 #gleaned from sam XT opt field
        self.stat_dict['n_of_adapter_clipped_reads'] = 0
        
        self.stat_dict['n_of_nt_covered'] = 0
        self.stat_dict['n_of_nt_termini'] = 0
        self.stat_dict['n_of_oneD_operations'] = 0
        
        return
    
    def __call__(self, ome_coords,
                       pysam_bam_file_conn,
                       cleaved_readid_db_conn = None,
                       all_adapter_clipped = False,
                       uniq_only = True,
                       oneD_rate_mode = True,
                       rate_cutoff = 15,
                       cleav_rate_mode = False,
                       stranded = True):
        
        reference, start, end, strand = ome_coords
        
        # fetch the bam-encoded alignment data for the region
        # delineated by the ome_coords.
        # note that it will also return reads that are only
        # partially overlapping with the region
        c = fetch_hitsclip_dict()
        pysam_bam_file_conn.fetch( reference,
                                   start,
                                   end,
                                   callback = c )
        fetched_alignments = c.l
        self.stat_dict['n_of_reads'] = len(fetched_alignments)
        
        # instantiate OmeDict, which will aggregate and then output the
        # relevant basewise data as it is gleaned from the alignment data
        raw_cover_cd = OmeDict()
        cleavage_cd = OmeDict()
        oneD_cd = OmeDict()
        
        for d in fetched_alignments:
            
            if strand == d['strand']:
                pass
            else:
                continue
            
            self.stat_dict['n_of_strand_matched_reads'] += 1
            
            #identify uniquely aligned reads by accessing the sam optdict
            if d['OPT_XT'] == 'U':
                self.stat_dict['n_of_unique_aligns'] += 1
            elif uniq_only:
                continue
            else:
                pass
            
            ####
            #store basewise raw coverage data in OmeDict
            raw_cover_st = SiteTuple(reference,
                                     start = d['POS'],
                                     end = d['AEND'],
                                     strand = d['strand'],
                                     data = 1)
            raw_cover_cd.add_site_tuple( raw_cover_st )
            
            ####
            #store basewise cleavage data in OmeDict (if present)
            self.adapter_clipped_bool = all_adapter_clipped
            if cleaved_readid_db_conn:
                if cleaved_readid_db_conn.readid_lookup(d['QNAME']):
                    self.adapter_clipped_bool = True
                    self.stat_dict['n_of_adapter_clipped_reads'] += 1
                else:
                    pass
            else:
                pass
            
            if self.adapter_clipped_bool or all_adapter_clipped:
                
                #add left cleavage site
                left_terminus_st = SiteTuple(reference,
                                             start = d['POS'],
                                             end = d['POS'] + 1,
                                             strand = d['strand'],
                                             data = 1)
                cleavage_cd.add_site_tuple( left_terminus_st )
                
                #add right cleavage site
                right_terminus_st = SiteTuple(reference,
                                              start = d['AEND'] - 1,
                                              end = d['AEND'],
                                              strand = d['strand'],
                                              data = 1)
                cleavage_cd.add_site_tuple( right_terminus_st )
            
            ####
            # store 1D (aka single nucleotide deletion) data in OmeDict
            # note: pysam uses an integer key to refer to CIGAR operations
            # http://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedRead.cigar
            # deletion operations are coded as a 2
            if (2, 1) in d['CIGAR']:
                # the postional index of the alignment
                # operation (on the reference) will be
                # stored as i
                i = -1
                for operation, length in d['CIGAR']:
                    if operation != 1:
                        i += length
                        if operation == 2 and length == 1:
                            oneD_st = SiteTuple(reference,
                                                start = d['POS'] + i,
                                                end = d['POS'] + i + 1,
                                                strand = d['strand'],
                                                data = 1)
                            oneD_cd.add_site_tuple(oneD_st)
                        else:
                            pass
        
        # output the data vectors. 
        # note that the get_data method reverses the bottom strand vectors by 
        # default (see orient_strands argument). Therefore the top and bottom 
        # strand coverage vectors have the same 5'->3' polarity.
        # also note that the stranded argument here can be toggled to accomodate
        # strand-agnostic sequencing libraries. However HITS-CLIP libraries are
        # typically stranded.
        if stranded:
            raw_cover_vec = raw_cover_cd.get_data(ome_coords)
            self.stat_dict['n_of_nt_covered'] = sum(raw_cover_vec)
            cleavage_vec = cleavage_cd.get_data(ome_coords)
            self.stat_dict['n_of_nt_termini'] = sum(cleavage_vec)
            oneD_vec = oneD_cd.get_data(ome_coords)
            self.stat_dict['n_of_oneD_operations'] = sum(oneD_vec)
        else:
            raw_cover_vec = raw_cover_cd.get_data(ome_coords, both_strands = True)
            self.stat_dict['n_of_nt_covered'] = sum(raw_cover_vec)
            cleavage_vec = cleavage_cd.get_data(ome_coords, both_strands = True)
            self.stat_dict['n_of_nt_termini'] = sum(cleavage_vec)
            oneD_vec = oneD_cd.get_data(ome_coords, both_strands = True)
            self.stat_dict['n_of_oneD_operations'] = sum(oneD_vec)
        
        if oneD_rate_mode == True:
            oneD_rate_vec = []
            for x, y in zip(oneD_vec, raw_cover_vec):
                if y < rate_cutoff:
                    oneD_rate_vec.append(0)
                else:
                    oneD_rate_vec.append(x/y)
            
            oneD_vec = oneD_rate_vec
        
        if cleav_rate_mode == True:
            cleav_rate_vec = []
            for x, y in zip(cleavage_vec, raw_cov_vec):
                if y < rate_cutoff:
                    cleav_rate_vec.append(0)
                else:
                    cleav_rate_vec.append(x/y)
            
            cleavage_vec = cleav_rate_vec
        
        vector_tuple = ( raw_cover_vec,
                         oneD_vec,
                         cleavage_vec )
        
        return vector_tuple

#NOTE: Called by ./hitsclip_bed-dump.py
class hitsclip_vectors_2_bg():
    '''
    Generate basewise counts for 1D, cleavage and raw coverage from the 
    set of alignments from the reference interval defined in the site_tup.
    This function parses the list of dictionaries generated by the pysam
    indexed bam file interface function at ??.
    Note: if the cleav_db_fp is set to None then it is assumed that all
    the reads in the sam file are adapter-clipped
    '''
    
    def __init__(self):
        
        self.stat_dict = {}
        self.stat_dict['n_of_mapped_reads'] = 0
        self.stat_dict['n_of_unique_aligns'] = 0 #gleaned from sam XT opt field
        self.stat_dict['n_of_adapter_clipped_reads'] = 0
        
        self.stat_dict['n_of_nt_covered'] = 0
        self.stat_dict['n_of_nt_termini'] = 0
        self.stat_dict['n_of_oneD_operations'] = 0
        
        return
    
    def __call__(self, pysam_bam_file_conn,
                       sample_name,
                       bg_fh_top_strand_cov,
                       bg_fh_bot_strand_cov,
                       bg_fh_top_strand_clv,
                       bg_fh_bot_strand_clv,
                       bg_fh_strand_1D,
                       bg_fh_strand_1D,
                       cleaved_readid_db_conn = None,
                       all_adapter_clipped = False,
                       uniq_only = True,
                       oneD_rate_mode = True,
                       oneD_rate_cutoff = 15,
                       cleav_rate_mode = False,
                       cleav_rate_cutoff = None,
                       max_chunk_size = 1000000):
        
        #see bedgraph definition at https://genome.ucsc.edu/goldenPath/help/bedgraph.html
        #note  bedGraphToBigWig can compress the file see http://genome.ucsc.edu/goldenPath/help/bigWig.html#Ex3
        generic_head = 'track type=bedGraph name="{0}" description="{1}{2}" visibility=full\n'
        bg_fh_top_strand_cov.write(generic_head.format('5\'--> coverage', sample_name, 'Top strand coverage'))
        bg_fh_bot_strand_cov.write(generic_head.format('<--3\' coverage', sample_name, 'Bottom strand coverage'))
        bg_fh_top_strand_term.write(generic_head.format('5\'--> termini', sample_name, 'Top strand termini'))
        bg_fh_bot_strand_term.write(generic_head.format('<--3\' coverage', sample_name, 'Bottom strand termini'))
        bg_fh_strand_1D.write(generic_head.format('5\'--> coverage', sample_name, 'Top strand single nucleotide deletions'))
        bg_fh_strand_1D.write(generic_head.format('<--3\' coverage', sample_name, 'Bottom single nucleotide deletions'))
        
        i = pysam_bam_file_conn.fetch()
        
        # instantiate OmeDict, which will aggregate and then output the
        # relevant basewise data as it is gleaned from the alignment data
        raw_cover_cd = OmeDict()
        cleavage_cd = OmeDict()
        oneD_cd = OmeDict()
        
        chunk_size = 0 # an object to manage memory tracks with n_of_nt_covered
        
        for i in fetched_alignments:
            
            self.stat_dict['n_of_mapped_reads'] += 1
            
            d = {}
            d['QNAME'] = alignment.qname
            d['POS'] = alignment.pos
            d['AEND'] = alignment.aend
            d['CIGAR'] = alignment.cigar
            d['IS_REVERSE'] = alignment.is_reverse
            if alignment.is_reverse:
                d['strand'] = '-'
            else:
                d['strand'] = '+'
            d['OPT_XT'] = alignment.opt('XT')
            
            #identify uniquely aligned reads by accessing the sam optdict
            if d['OPT_XT'] == 'U':
                self.stat_dict['n_of_unique_aligns'] += 1
            elif uniq_only:
                continue
            else:
                pass
            
            ####
            #store basewise raw coverage data in OmeDict
            raw_cover_st = SiteTuple(reference,
                                     start = d['POS'],
                                     end = d['AEND'],
                                     strand = d['strand'],
                                     data = 1)
            n = raw_cover_cd.add_site_tuple( raw_cover_st )
            self.stat_dict['n_of_nt_covered'] += n
            chunk_size += n
            
            ####
            #store basewise cleavage data in OmeDict (if present)
            self.adapter_clipped_bool = all_adapter_clipped
            if cleaved_readid_db_conn:
                if cleaved_readid_db_conn.readid_lookup(d['QNAME']):
                    self.adapter_clipped_bool = True
                    self.stat_dict['n_of_adapter_clipped_reads'] += 1
                else:
                    pass
            else:
                pass
            
            if self.adapter_clipped_bool or all_adapter_clipped:
                
                #add left cleavage site
                left_terminus_st = SiteTuple(reference,
                                             start = d['POS'],
                                             end = d['POS'] + 1,
                                             strand = d['strand'],
                                             data = 1)
                n = cleavage_cd.add_site_tuple( left_terminus_st )
                
                #add right cleavage site
                right_terminus_st = SiteTuple(reference,
                                              start = d['AEND'] - 1,
                                              end = d['AEND'],
                                              strand = d['strand'],
                                              data = 1)
                n = cleavage_cd.add_site_tuple( right_terminus_st )
            
            ####
            # store 1D (aka single nucleotide deletion) data in OmeDict
            # note: pysam uses an integer key to refer to CIGAR operations
            # http://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedRead.cigar
            # deletion operations are coded as a 2
            if (2, 1) in d['CIGAR']:
                # the postional index of the alignment
                # operation (on the reference) will be
                # stored as i
                i = -1
                for operation, length in d['CIGAR']:
                    if operation != 1:
                        i += length
                        if operation == 2 and length == 1:
                            oneD_st = SiteTuple(reference,
                                                start = d['POS'] + i,
                                                end = d['POS'] + i + 1,
                                                strand = d['strand'],
                                                data = 1)
                            n = oneD_cd.add_site_tuple(oneD_st)
                        else:
                            pass
            
            if chunk_size >= max_chunk_size:
                #TODO:writeout all coords 5' of the most recent alignment and delete
                # the key-value pairs from memory. include a print statement that 
                # declares the writeout is occuring so that memory clean-up
                # can be verified by watching htop
                
                # output the data vectors. 
                # note that the get_data method reverses the bottom strand vectors by 
                # default (see orient_strands argument). Therefore the top and bottom 
                # strand coverage vectors have the same 5'->3' polarity.
                
                # NOTE: the coordinates of the current alignment comprises
                # the boundary coords. Eveerything upstream will be dumped
                # to file and deleted from memory.
                
                n = cleavage_cd.dump2bg( bg_fh_top_strand_term, 
                                         bg_fh_bot_strand_term, 
                                         ome_coords, # current boundary coords
                                         rate_mode = cleav_rate_mode, 
                                         rate_cutoff = cleav_rate_cutoff,
                                         rate_denom_cd = raw_cover_cd,
                                         purge_memory = True,
                                        )
                self.stat_dict['n_of_nt_termini'] += n
                
                n = oneD_cd.dump2bg( bg_fh_top_strand_1D, 
                                     bg_fh_bot_strand_1D, 
                                     ome_coords, # current boundary coords
                                     rate_mode = oneD_rate_mode, 
                                     rate_cutoff = oneD_rate_cutoff,
                                     rate_denom_cd = raw_cover_cd,
                                     purge_memory = True,
                                    )
                self.stat_dict['n_of_oneD_operations'] += n
                
                # The raw coverage ome_dict must be dumped last because it 
                # serves as the denominator for the rate calculations for 
                # cleavage and oneD
                n = raw_cover_cd.dump2bg( bg_fh_top_strand_cov, 
                                          bg_fh_bot_strand_cov, 
                                          ome_coords, # current boundary coords
                                          purge_memory = True,
                                        )
                self.stat_dict['n_of_nt_covered'] += n
        
        return

#TODO: implement me
class build_parclip_vectors():
    '''
    Generate basewise counts for 1D, cleavage and raw coverage from the 
    set of alignments from the reference interval defined in the site_tup.
    This function parses the list of dictionaries generated by the pysam
    indexed bam file interface function at ??.
    Note: if the cleav_db_fp is set to None then it is assumed that all
    the reads in the sam file are adapter-clipped
    '''
    
    def __init__(self):
        
        self.stat_dict = {}
        self.stat_dict['n_of_reads'] = 0
        self.stat_dict['n_of_strand_matched_reads'] = 0
        self.stat_dict['n_of_unique_aligns'] = 0 #gleaned from sam XT opt field
        self.stat_dict['n_of_adapter_clipped_reads'] = 0
        
        self.stat_dict['n_of_nt_covered'] = 0
        self.stat_dict['n_of_nt_termini'] = 0
        self.stat_dict['n_of_oneD_operations'] = 0
        
        return
    
    def __call__(self, ome_coords,
                       pysam_bam_file_conn,
                       cleaved_readid_db_conn = None,
                       all_adapter_clipped = False,
                       uniq_only = True,
                       oneD_rate_mode = True,
                       rate_cutoff = 15,
                       cleav_rate_mode = False,
                       stranded = True):
        
        reference, start, end, strand = ome_coords
        
        # fetch the bam-encoded alignment data for the region
        # delineated by the ome_coords.
        # note that it will also return reads that are only
        # partially overlapping with the region
        c = fetch_hitsclip_dict()
        pysam_bam_file_conn.fetch( reference,
                                   start,
                                   end,
                                   callback = c )
        fetched_alignments = c.l
        self.stat_dict['n_of_reads'] = len(fetched_alignments)
        
        # instantiate OmeDict, which will aggregate and then output the
        # relevant basewise data as it is gleaned from the alignment data
        raw_cover_cd = OmeDict()
        cleavage_cd = OmeDict()
        oneD_cd = OmeDict()
        
        for d in fetched_alignments:
            
            if strand == d['strand']:
                pass
            else:
                continue
            self.stat_dict['n_of_strand_matched_reads'] += 1
            
            #identify uniquely aligned reads by accessing the sam optdict
            if d['OPT_XT'] == 'U':
                self.stat_dict['n_of_unique_aligns'] += 1
            elif uniq_only:
                continue
            else:
                pass
            
            ####
            #store basewise raw coverage data in OmeDict
            raw_cover_st = SiteTuple(reference,
                                     start = d['POS'],
                                     end = d['AEND'],
                                     strand = d['strand'],
                                     data = 1)
            raw_cover_cd.add_site_tuple( raw_cover_st )
            
            ####
            #store basewise cleavage data in OmeDict (if present)
            self.adapter_clipped_bool = all_adapter_clipped
            if cleaved_readid_db_conn:
                if cleaved_readid_db_conn.readid_lookup(d['QNAME']):
                    self.adapter_clipped_bool = True
                    self.stat_dict['n_of_adapter_clipped_reads'] += 1
                else:
                    pass
            else:
                pass
            
            if self.adapter_clipped_bool or all_adapter_clipped:
                
                #add left cleavage site
                left_terminus_st = SiteTuple(reference,
                                             start = d['POS'],
                                             end = d['POS'] + 1,
                                             strand = d['strand'],
                                             data = 1)
                cleavage_cd.add_site_tuple( left_terminus_st )
                
                #add right cleavage site
                right_terminus_st = SiteTuple(reference,
                                              start = d['AEND'] - 1,
                                              end = d['AEND'],
                                              strand = d['strand'],
                                              data = 1)
                cleavage_cd.add_site_tuple( right_terminus_st )
            
            ####
            # store 1D (aka single nucleotide deletion) data in OmeDict
            # note: pysam uses an integer key to refer to CIGAR operations
            # http://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedRead.cigar
            # deletion operations are coded as a 2
            if (2, 1) in d['CIGAR']:
                # the postional index of the alignment
                # operation (on the reference) will be
                # stored as i
                i = -1
                for operation, length in d['CIGAR']:
                    if operation != 1:
                        i += length
                        if operation == 2 and length == 1:
                            oneD_st = SiteTuple(reference,
                                                start = d['POS'] + i,
                                                end = d['POS'] + i + 1,
                                                strand = d['strand'],
                                                data = 1)
                            oneD_cd.add_site_tuple(oneD_st)
                        else:
                            pass

#def detect_mm_variant(self, native_nt = 'T', variant_nt = 'C', 
                            #siteTup_mode = False, 
                            #nt_set = 'ATCGN'):
    #'''
    #you can check a specific nt variant.
    #the default mismatch that is detected is T to C.
    #returns a list of integers for each mismatch event.
    #the integers returned are the zero-based indices of the variant 
    #occurence along the string stored in SEQ field.
    #a list of zero-based chromosomal coordinates can be returned instead 
    #by setting siteTup_mode to True. the chromosomal coordinates will
    #be returned as a list of siteTups can be calculated.
    #'''
    #v = []
    ##NOTE: this read must be mapped. otherwise optdict will be empty.
    #sam_data_type, mismatch_str = self.d['OPTDICT']['MD']
    #if native_nt in mismatch_str:
        #m = self.MDZregexP.findall(mismatch_str)
        #sam_seq_pos = -1 # this is the zero-based index into the SEQ string
        #ref_pos = self.d['POS'] - 2 # this creates the zero-based index into the REF
        #for i in m:
            #if native_nt in i and '^' != i[0]:
                #for nt in i:
                    #sam_seq_pos += 1
                    #ref_pos += 1
                    #if nt == native_nt:
                        #if self.d['SEQ'][sam_seq_pos] == variant_nt:
                            ##print('found one', self.d['QNAME']) #debugging
                            #if siteTup_mode:
                                #v.append( (self.d['RNAME'], ref_pos, ref_pos + 1, self.gleanStrand()) )
                            #else:
                                #v.append( (sam_seq_pos, ref_pos) )
                        #else:
                            #continue
                    #else:
                        #continue
            
            #elif '^' == i[0]:
                #ref_pos -= len(i[1:])
                #continue
            
            #elif i in nt_set:
                #sam_seq_pos += len(i)
                #ref_pos += len(i)
                #continue
            
            #elif re.match('[0-9]+', i):
                #sam_seq_pos += int(i)
                #ref_pos += int(i)
                #continue
            
            #else:
                #print('something went wrong.')
                #raise IOError
    
    #return v


        ## output the data vectors. 
        ## note that the get_data method reverses the bottom strand vectors by 
        ## default (see orient_strands argument). Therefore the top and bottom 
        ## strand coverage vectors have the same 5'->3' polarity.
        ## also note that the stranded argument here can be toggled to accomodate
        ## strand-agnostic sequencing libraries. However HITS-CLIP libraries are
        ## typically stranded.
        #if stranded:
            #raw_cover_vec = raw_cover_cd.get_data(ome_coords)
            #self.stat_dict['n_of_nt_covered'] = sum(raw_cover_vec)
            #cleavage_vec = cleavage_cd.get_data(ome_coords)
            #self.stat_dict['n_of_nt_termini'] = sum(cleavage_vec)
            #oneD_vec = oneD_cd.get_data(ome_coords)
            #self.stat_dict['n_of_oneD_operations'] = sum(oneD_vec)
        #else:
            #raw_cover_vec = raw_cover_cd.get_data(ome_coords, both_strands = True)
            #self.stat_dict['n_of_nt_covered'] = sum(raw_cover_vec)
            #cleavage_vec = cleavage_cd.get_data(ome_coords, both_strands = True)
            #self.stat_dict['n_of_nt_termini'] = sum(cleavage_vec)
            #oneD_vec = oneD_cd.get_data(ome_coords, both_strands = True)
            #self.stat_dict['n_of_oneD_operations'] = sum(oneD_vec)
        
        #if oneD_rate_mode == True:
            #oneD_rate_vec = []
            #for x, y in zip(oneD_vec, raw_cover_vec):
                #if y < rate_cutoff:
                    #oneD_rate_vec.append(0)
                #else:
                    #oneD_rate_vec.append(x/y)
            
            #oneD_vec = oneD_rate_vec
        
        #if cleav_rate_mode == True:
            #cleav_rate_vec = []
            #for x, y in zip(cleavage_vec, raw_cov_vec):
                #if y < rate_cutoff:
                    #cleav_rate_vec.append(0)
                #else:
                    #cleav_rate_vec.append(x/y)
            
            #cleavage_vec = cleav_rate_vec
        
        #vector_tuple = ( raw_cover_vec,
                         #oneD_vec,
                         #cleavage_vec )
        
        #return vector_tuple

##TESTING
##http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr11%3A118964585-118966177&hgsid=369409831_ZPHI1ErsOPVE9y4sB1fWqubbTvzV
##note these are hg19 coords
#ome_coords = ('chr11', 118964584, 118966177, '-')
#bam_fp = '/home/lbthrice/Desktop/data_sample/Brooks_SLBP/s_1xS01_sequence.PP.bwa-hg19.bam'
#cleaved_readid_db_fp = 's_1xS01_sequence.PP.readids'
#with pysam.Samfile( bam_fp, "rb" ) as bam_conn, \
#     ReadidSQLite( cleaved_readid_db_fp ) as readid_conn:
#    
#    ome_coords = ref, start, end, strand
#    
#    f = build_hitsclip_vectors()
#    t = f( ome_coords,
#           bam_conn, 
#           cleaved_readid_db_conn = readid_conn )
#    
#    raw_cover_l, cleavage_l, oned_rate_l = t
#    
#    f.stat_dict
##Out[2]: 
##{'n_of_adapter_clipped_reads': 431,
## 'n_of_nt_covered': 24749,
## 'n_of_nt_termini': 862,
## 'n_of_oneD_operations': 45,
## 'n_of_reads': 1048,
## 'n_of_strand_matched_reads': 1047,
## 'n_of_unique_aligns': 895}
