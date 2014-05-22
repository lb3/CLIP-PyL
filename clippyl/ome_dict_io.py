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
        """
        reference, start, end, strand, data = site_tuple
        
        if reference in self.d[strand]:
            for coord in range(start, end):
                if coord in self.d[strand][reference]:
                    self.d[strand][reference][coord] += data
                else:
                    self.d[strand][reference][coord] = data
        else:
            self.d[strand][reference] = OrderedDict()
            for coord in range(start, end):
                self.d[strand][reference][coord] = data
        
        return
    
    def get_data(self, ome_coords, orient_strands = True):
        """
        Retrieve the data corresponding to the given coordinates from 
        the OmeDict. A list of integers is returned (aka a coverage vector).
        The orient_strands argument caused the minus strand data lists to
        reverse such that the top and bottom strand vectors are oriented in the
        same left to right (5'->3') polarity.
        """
        reference, start, end, strand = ome_coords
        
        data = []
        for coord in range(start, end):
            data.append(self.d[strand][reference].get(coord, 0))
        
        if orient_strands and strand == '-':
            data.reverse()
        
        return data
