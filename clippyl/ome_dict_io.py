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
    def dump2bg(self, bg_fh_top_strand, 
                      bg_fh_bot_strand, 
                      b_ref, 
                      b_coord, 
                      rate_mode = False, # use rate mode for normalization
                      rate_cutoff = 0, # minimum number of sequenced bases required to calculate a coverage rate
                      rate_denom_od = {}, # ome dict containing rate denominator (the numerator is the ome_dict)
                      ):
        
        # the function will return these coverage stats
        n_of_nt_covered = 0
        n_of_sequenced_bases = 0
        
        # note: I used nested functions for writing. This allowed me to use 
        # return to break out of the nested loop
        # http://stackoverflow.com/a/189685
        # https://mail.python.org/pipermail/python-3000/2007-July/008663.html
        def write_strand(b_ref, b_coord, strand, bg_fh):
            # write top strand until the boundary coords are encountered
            # note that dictionary iteration methods have distinct semantics
            # previous python versions
            # http://stackoverflow.com/a/6777632/892030
            # http://legacy.python.org/dev/peps/pep-0469/
            
            for ref in list(self.d[strand].keys()):
                
                coord_l = list(self.d[strand][ref])
                
                # get the first coord value for this reference
                i = 0
                coord = coord_l[i]
                datum = self.d[strand][ref][coord_l[i]]
                assert type(coord) == int
                assert type(datum) == int
                # initial check that we are not beyond the boundary
                if b_ref == ref and b_coord <= coord:
                    # we are beyond the boundary, exit
                    return
                else:
                    # begin dumping bases in this reference...
                    # first instantiate objects to hold values for 
                    # write out to bedgraph
                    ref_writeout = ref
                    start_coord_writeout = coord
                    end_coord_writeout = coord + 1
                    n_of_nt_covered += 1
                    n_of_sequenced_bases += datum
                    
                    if rate_mode:
                        # with rate_cutoff nothing will be written if there are
                        # not enough sequenced bases in the denominator
                        # otherwise, it is assumed that there are enough 
                        # observations to calculate a rate
                        if rate_denom_od[strand][ref][coord] < rate_cutoff
                            datum_writeout = None
                        else:
                            datum_writeout = datum/rate_denom_od[strand][ref][coord]
                            # store the previous values to allow sequential coords that
                            # hold the same value to be written on a single line (as a 
                            # range) in accordance with the bedgraph format
                            prev_coord = coord
                            prev_datum = datum/rate_denom_od[strand][ref][coord]
                    else:
                        datum_writeout = datum
                        # store the previous values to allow sequential coords that
                        # hold the same value to be written on a single line (as a 
                        # range) in accordance with the bedgraph format
                        prev_coord = coord
                        prev_datum = datum
                    
                    # purge from memory
                    del self.d[strand][ref][coord]
                    
                    if len(coord_l) > 1:
                        for coord in coord_l[1:]:
                            if b_ref == ref and b_coord <= coord:
                                # past boundary, write current values to file
                                # and exit func
                                if datum_writeout != None:
                                    bg_fh.write('\t'.join(ref_writeout,
                                                          str(start_coord_writeout),
                                                          str(end_coord_writeout),
                                                          '{:f}'.format(datum_writeout)) + '\n')
                                    return n_of_sequenced_bases, n_of_nt_covered
                                else:
                                    return n_of_sequenced_bases, n_of_nt_covered
                            else:
                                pass
                            
                            datum = self.d[strand][ref][coord]
                            del self.d[strand][ref][coord]
                            n_of_nt_covered += 1
                            n_of_sequenced_bases += datum
                            
                            if rate_mode:
                                # with rate_cutoff nothing will be written if there are
                                # not enough sequenced bases in the denominator
                                # otherwise, it is assumed that there are enough 
                                # observations to calculate a rate
                                if rate_denom_od[strand][ref][coord] < rate_cutoff
                                    datum = None
                                else:
                                    datum = datum/rate_denom_od[strand][ref][coord]
                            else:
                                pass
                            
                            # convert to string representation to avoid floating
                            # point arithmetic issues
                            # https://docs.python.org/3.3/tutorial/floatingpoint.html#floating-point-arithmetic-issues-and-limitations
                            a = '{:f}'.format(datum)
                            b = '{:f}'.format(prev_datum)
                            if coord == prev_coord + 1 and  a == b:
                                # combine sequential equivalent bases for writeout 
                                # as a range in the bedgraph format
                                # update placeholder objects
                                end_coord_writeout = coord + 1
                                prev_coord = coord
                                
                            else:
                                
                                bg_fh.write('\t'.join(ref_writeout = ref,
                                                      start_coord_writeout,
                                                      end_coord_writeout,
                                                      '{:f}'.format(datum_writeout)) + '\n')
                                
                                # update objects to hold values for write out to bedgraph
                                start_coord_writeout = coord
                                end_coord_writeout = coord + 1
                                datum_writeout = datum
                                
                                # store the previous values to allow sequential coords that
                                # hold the same value to be written on a single line (as a 
                                # range) in the bedgraph format
                                prev_coord = coord
                                prev_datum = datum
                                
                    
                    # final writout for this reference
                    bg_fh.write('\t'.join(ref_writeout = ref,
                                          start_coord_writeout,
                                          end_coord_writeout,
                                          {':f'}.format(datum_writeout)) + '\n')
                    
                    return n_of_nt_covered, n_of_sequenced_bases
        
        t = write_strand(b_ref = b_ref, 
                         b_coord = b_coord, 
                         strand = '+', 
                         bg_fh = bg_fh_top_strand)
        
        print(t) # debugging
        n_of_nt_covered, n_of_sequenced_bases = t
        
        t = write_strand(b_ref = b_ref, 
                         b_coord = b_coord, 
                         strand = '-', 
                         bg_fh = bg_fh_bot_strand)
        
        n_of_nt_covered += t[0]
        n_of_sequenced_bases += t[1]
        print(t) # debugging
        
        return n_of_nt_covered, n_of_sequenced_bases

