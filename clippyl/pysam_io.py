import pysam
#docs:
#pysam development branch: https://github.com/pysam-developers/pysam
#pysam docs: http://pysam.readthedocs.org/en/latest/
#pysam docs hosted here too: http://www.cgat.org/~andreas/documentation/pysam/api.html

class fetch_hitsclip_dict:
    
    '''This is a callback function for use with the pysam fetch method. It
       fetches the alignment information necessary to build the hitsclip
       signature.'''
    
    def __init__(self):
        self.l = []
        return
    
    def __enter__(self):
        return self
    
    def __exit__(self, *exc_info):
        self.fh.close()
    
    def __call__(self, alignment):
        #print(type(alignment))
        ##NOTE: the alignment object is <class 'pysam.csamtools.AlignedRead'>
        ##see API documentation for a list if attributes: 
        ##https://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedRead
        
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
        self.l.append(d)
        
        return

class fetch_iclip_dict:
    
    '''This is a callback function for use with the pysam fetch method. It
       fetches the alignment information necessary to build the iclip
       signature.'''
    
    def __init__(self):
        self.l = []
        return
    
    def __enter__(self):
        return self
    
    def __exit__(self, *exc_info):
        self.fh.close()
    
    def __call__(self, alignment):
        #print(type(alignment))
        ##NOTE: the alignment object is <class 'pysam.csamtools.AlignedRead'>
        ##see API documentation for a list if attributes: 
        ##https://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedRead
        
        d = {}
        d['QNAME'] = alignment.qname # Name of read
        d['POS'] = alignment.pos
        d['AEND'] = alignment.aend
        d['IS_REVERSE'] = alignment.is_reverse
        if alignment.is_reverse:
            d['strand'] = '-'
        else:
            d['strand'] = '+'
        d['OPT_XT'] = alignment.opt('XT')
        self.l.append(d)
        
        return

class fetch_parclip_dict:
    
    '''This is a callback function for use with the pysam fetch method. It
       fetches the alignment information necessary to build the parclip
       signature.'''
    
    def __init__(self):
        self.l = []
        return
    
    def __enter__(self):
        return self
    
    def __exit__(self, *exc_info):
        self.fh.close()
    
    def __call__(self, alignment):
        #print(type(alignment))
        ##NOTE: the alignment object is <class 'pysam.csamtools.AlignedRead'>
        ##see API documentation for a list if attributes: 
        ##https://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedRead
        
        d = {}
        d['QNAME'] = alignment.qname # Name of read
        d['POS'] = alignment.pos
        d['AEND'] = alignment.aend
        d['IS_REVERSE'] = alignment.is_reverse
        if alignment.is_reverse:
            d['strand'] = '-'
        else:
            d['strand'] = '+'
        d['IS_REVERSE'] = alignment.is_reverse
        d['OPT_XT'] = alignment.opt('XT')
        d['OPT_MD'] = alignment.opt('MD')
        self.l.append(d)
        
        return

