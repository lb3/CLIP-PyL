from collections import namedtuple, OrderedDict

#Custom namedtuple class for ChromDict transactions
SiteTuple = namedtuple('SiteTuple', ['reference', 
                                     'start', 
                                     'end', 
                                     'strand', 
                                     'data'])

class OmeDict():
    """
    An dictionary-based storage structure for "omics" data (hence "OmeDict").
    The supported input object type is the SiteTuple, which is a namedtuple
    with the following fields:
    reference, start, end, strand, dataobject
    
    The reference field is most often populated by a chromosome name but, 
    alternatively, a transcript name or some other reference type is also valid.
    The start and end fields are coordinates along the reference. The 
    coordinates must be "UCSC-style" (right-open zero based; (start, end]). Note
    that the supported object type for the data field is integer and that the 
    addition operator is used to add data.
    
    Add data to the OmeDict using the add_site_tuple method.
    
    Also, note that an OrderedDict is used to store the references as well as
    the coordinate keys for each reference. This makes retrieving data from the
    OmeDict more convinient; key input order is preserved, which obviates the
    need to sort downstream. WARNING: This function is appropriate when the 
    input is pre-sorted (e.g. sam and bed files).
    """
    
    def __init__(self):
        self.d = {'+' : OrderedDict(), '-' : OrderedDict()}
        return
    
    def add_site_tuple(self, site_tuple):
        """
        The supported input object type is the SiteTuple class, which is a 
        namedtuple with the following fields: reference, start, end, strand, 
        dataobject. Note that the supported object type for the data
        field is integer and that the addition operator is used to add data.
        The number of new coords written is returned.
        """
        reference, start, end, strand, data = site_tuple
        
        n_new_coord_keys = 0
        
        if reference in self.d[strand]:
            for coord in range(start, end):
                if coord in self.d[strand][reference]:
                    self.d[strand][reference][coord] += data
                else:
                    self.d[strand][reference][coord] = data
                    n_new_coord_keys += 1
        else:
            self.d[strand][reference] = OrderedDict()
            for coord in range(start, end):
                self.d[strand][reference][coord] = data
                n_new_coord_keys += 1
        
        return n_new_coord_keys
    
    def get_data(self, ome_coords, both_strands = False, orient_strands = True):
        """
        Retrieve the data corresponding to the given coordinates from 
        the OmeDict. A list of integers is returned (aka a coverage vector).
        The orient_strands argument caused the minus strand data lists to
        reverse such that the top and bottom strand vectors are oriented in the
        same left to right (5'->3') polarity. The both_strands argument here 
        can be toggled to accomodate strand-agnostic sequencing libraries. 
        However HITS-CLIP libraries are typically stranded.
        """
        reference, start, end, strand = ome_coords
        
        if both_strands:
            t_strand_data = []
            b_strand_data = []
            for coord in range(start, end):
                try:
                    t_strand_data.append(self.d['+'][reference].get(coord, 0))
                except KeyError:
                    t_strand_data.append(0)
                try:
                    b_strand_data.append(self.d['-'][reference].get(coord, 0))
                except KeyError:
                    b_strand_data.append(0)
            data = [sum(x) for x in zip(t_strand_data, b_strand_data)]
        else:
            data = []
            for coord in range(start, end):
                try:
                    data.append(self.d[strand][reference].get(coord, 0))
                except KeyError:
                    data.append(0)
        
        if orient_strands and strand == '-':
            data.reverse()
        
        return data

#TODO: improve efficiency by deleting dict in single operation
def dump2bg(ome_dict,
              bg_fh_top_strand, 
              bg_fh_bot_strand, 
              b_ref, 
              b_coord, 
              rate_mode = False, # use rate mode for normalization
              rate_cutoff = 0, # minimum number of sequenced bases required to calculate a coverage rate
              rate_denom_cd = {}, # ome dict containing rate denominator (the numerator is the ome_dict)
              ):
    
    # note: I used nested functions for writing. This allowed me to use 
    # return to break out of the nested loop
    # http://stackoverflow.com/a/189685
    # https://mail.python.org/pipermail/python-3000/2007-July/008663.html
    class write_strand():
        
        def __init__(self):
            
            # write top strand until the boundary coords are encountered
            # note that dictionary iteration methods have distinct semantics
            # previous python versions
            # http://stackoverflow.com/a/6777632/892030
            # http://legacy.python.org/dev/peps/pep-0469/
            
            # the function will return these coverage stats
            self.n_of_nt_covered = 0
            self.n_of_sequenced_bases = 0
            
            self.data_writeout = None
            self.ref_writeout =None
            self.start_coord_writeout = None
            self.end_coord_writeout = None
            self.datum_writeout = None
        
        def __call__(self, b_ref, b_coord, strand, bg_fh):
            
            coord = None
            datum = None
            
            for ref in list(ome_dict.d[strand].keys()):
                
                print('first ref', ref)
                
                coord_l = list(ome_dict.d[strand][ref].keys())
                print(len(coord_l), 'coords found')
                if not coord_l:
                    print('empty ref', ref) #debugging
                    del ome_dict.d[strand][ref]
                    continue
                coord = coord_l[0]
                
                # get the first coord value for this reference
                datum = ome_dict.d[strand][ref][coord]
                assert type(coord) == int
                print('first coord', coord)
                assert type(datum) == int
                print('first datum', datum)
                assert type(b_coord == int)
                # initial check that we are not beyond the boundary
                if b_ref == ref and b_coord <= coord:
                    # we are beyond the boundary, exit
                    print('initial boundary ref', ref) # debugging
                    print('initial boundary coord', coord) # debugging
                    return self.n_of_sequenced_bases, self.n_of_nt_covered
                else:
                    # begin dumping bases in this reference...
                    # first instantiate objects to hold values for 
                    # write out to bedgraph
                    # store the previous values to allow sequential coords that
                    # hold the same value to be written on a single line (as a 
                    # range) in accordance with the bedgraph format
                    self.ref_writeout = ref
                    self.start_coord_writeout = coord
                    self.end_coord_writeout = coord + 1
                    
                    self.n_of_nt_covered += 1
                    self.n_of_sequenced_bases += datum
                    
                    #prev_coord = coord
                    
                    if rate_mode:
                        # with rate_cutoff nothing will be written if there are
                        # not enough sequenced bases in the denominator
                        # otherwise, it is assumed that there are enough 
                        # observations to calculate a rate
                        if rate_denom_cd.d[strand][ref][coord] < rate_cutoff:
                            self.datum_writeout = None
                            # store the previous values to allow sequential coords that
                            # hold the same value to be written on a single line (as a 
                            # range) in accordance with the bedgraph format
                            #prev_datum = None
                        else:
                            self.datum_writeout = datum/rate_denom_cd.d[strand][ref][coord]
                            # store the previous values to allow sequential coords that
                            # hold the same value to be written on a single line (as a 
                            # range) in accordance with the bedgraph format
                            #prev_datum = datum/rate_denom_cd.d[strand][ref][coord]
                    else:
                        self.datum_writeout = datum
                        # store the previous values to allow sequential coords that
                        # hold the same value to be written on a single line (as a 
                        # range) in accordance with the bedgraph format
                        #prev_datum = datum
                    
                    # purge from memory
                    del ome_dict.d[strand][ref][coord]
                    
                    for coord in coord_l[1:]:
                        
                        if b_ref == ref and b_coord <= coord:
                            # past boundary, write current values to file
                            # and exit func
                            print('boundary ref', ref) # debugging
                            print('boundary coord', coord) # debugging
                            if self.datum_writeout != None:
                                bg_fh.write('\t'.join([self.ref_writeout,
                                                       str(self.start_coord_writeout),
                                                       str(self.end_coord_writeout),
                                                       '{:f}'.format(self.datum_writeout)]) + '\n')
                                return self.n_of_sequenced_bases, self.n_of_nt_covered
                            else:
                                return self.n_of_sequenced_bases, self.n_of_nt_covered
                        else:
                            pass
                        
                        datum = ome_dict.d[strand][ref][coord]
                        del ome_dict.d[strand][ref][coord]
                        self.n_of_nt_covered += 1
                        self.n_of_sequenced_bases += datum
                        
                        #TODO: make this into a method, tis redundant with above
                        if rate_mode:
                            # with rate_cutoff nothing will be written if there are
                            # not enough sequenced bases in the denominator
                            # otherwise, it is assumed that there are enough 
                            # observations to calculate a rate
                            if rate_denom_cd.d[strand][ref][coord] < rate_cutoff:
                                datum = None
                            else:
                                datum = datum/rate_denom_cd.d[strand][ref][coord]
                        else:
                            pass
                        
                        # convert to string representation to avoid floating
                        # point arithmetic issues
                        # https://docs.python.org/3.3/tutorial/floatingpoint.html#floating-point-arithmetic-issues-and-limitations
                        #print('datum', datum)
                        #print('previous datum', prev_datum)
                        
                        if datum != None and self.datum_writeout != None:
                            a = '{:f}'.format(datum)
                            b = '{:f}'.format(self.datum_writeout)
                            if coord == self.end_coord_writeout and a == b:
                                # combine sequential equivalent bases for writeout 
                                # as a range in the bedgraph format
                                # update placeholder objects
                                self.end_coord_writeout = coord + 1
                                #prev_coord = coord
                                continue
                        
                        if self.datum_writeout != None:
                            
                            bg_fh.write('\t'.join([self.ref_writeout,
                                                   str(self.start_coord_writeout),
                                                   str(self.end_coord_writeout),
                                                   '{:f}'.format(self.datum_writeout)]) + '\n')
                        
                        # update objects to hold values for write out to bedgraph
                        self.start_coord_writeout = coord
                        self.end_coord_writeout = coord + 1
                        self.datum_writeout = datum
                        
                        # store the previous values to allow sequential coords that
                        # hold the same value to be written on a single line (as a 
                        # range) in the bedgraph format
                        #prev_coord = coord
                        #prev_datum = datum
                        
                        
                    # final writout for this reference
                    if self.data_writeout != None:
                        
                        bg_fh.write('\t'.join([self.ref_writeout,
                                               str(self.start_coord_writeout),
                                               str(self.end_coord_writeout),
                                               '{:f}'.format(self.datum_writeout)]) + '\n')
            
            print('done')
            return self.n_of_nt_covered, self.n_of_sequenced_bases
    
    print('dumping top strand') #debugging
    f = write_strand()
    t = f(b_ref, 
          b_coord, 
          '+', 
          bg_fh_top_strand)
    print('dumped top strand',t) # debugging
    n_of_nt_covered, n_of_sequenced_bases = t
    
    print('dumping bottom strand') #debugging
    f = write_strand()
    t = f(b_ref, 
          b_coord, 
          '-', 
          bg_fh_bot_strand)
    print('dumped bottom strand',t) # debugging
    n_of_nt_covered += t[0]
    n_of_sequenced_bases += t[1]
    
    return n_of_nt_covered, n_of_sequenced_bases

