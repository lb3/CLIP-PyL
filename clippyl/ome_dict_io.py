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
    need to sort downstream.
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
    
    #TODO: build me purge memory in a single operation by reinitiating the 
    #      ome_dict at self and writing the data from beyong the boundary coords
    
    #TODO: use get_data??
    #NOTE: called by .vector_factory.hitsclip_vectors_2_bg
    def dump2bg(self, bg_fh_top_strand_1D, 
                      bg_fh_bot_strand_1D, 
                      ome_coords, 
                      rate_mode = False, 
                      rate_cutoff = None, 
                      rate_denom_cd = None, # ome dict containing rate denominator is typically raw coverage
                      purge_memory = True ):
        
        sum(raw_cover_vec) #TODO: return n_of_nt_covered
        
        # write top strand
            for ref in self.d['+']:
                coord, datum = self.d['+'][ref].items()[0]
                while coord, datum in self.d['+'][ref].items()[1:]:
                    
                    datum = self.d['+'][ref][coord]
                    curr_site_tup = (ref, coord, datum)
                    if coord ==
                    if r == boundary_ref and c > boundary_coord:
                       
                        return
                    else:
                        bg_fh.write('\t'.join() + '\n')
        
        # write bottom strand
            for r in self.d['-']:
                for c in self.d['-'][r]:
                   if r == boundary_ref and c == boundary_coord:
                       return
        
        reference, start, end, strand, data = site_tuple
            
ome_coords, both_strands = False, orient_strands = True):
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

