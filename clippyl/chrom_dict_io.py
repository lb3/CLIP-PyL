

#Custom namedtuple class for ChromDict transactions
SiteTuple = namedtuple('SiteTuple', ['reference', 
                                     'start', 
                                     'end', 
                                     'strand', 
                                     'data'])

class ChromDict():
    """
    An dictionary-based storage structure for storing data about
    chromosome positions. The supported input object type is the SiteTuple
    class, which is a namedtuple with the following fields:
    (str(reference), int(startCoord), int(endCoord), str(strand), dataobject)
    
    The reference field is most often populated by a chromosome name but, 
    alternatively, a gene model name or some other reference is also valid. The
    coordinates expected are UCSC-style (right-open zero based; (start, end]).
    
    Add data to the ChromDict using the add_site_tuple method. Note that the
    supported object type for the data field is integer and that the addition
    operator is used to add data. Also, note that an OrderedDict is used to
    store the references as well as the coordinate keys for each reference.
    This perserved input order and makes retrieving the data easier, obviating
    the need to sort downstream.
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

