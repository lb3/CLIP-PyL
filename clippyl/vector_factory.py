

#this function takes a tup of clipPyL db that hav been loaded
#NOTE: template code is graphic_mpl.iCliPyL_plot
#NOTE: negative strand vectors are autromatically reversed 
# so that all outout graphs are in the same orientation.
#TODO: parameterize negative strand reversal upstream at 
# chrom_db_io.siteTup2sigTup
#TODO: change cleavCov object name for iCLIP types; it is not cleav it is RT term
def iCLIPvectors(   bed_gen, 
                    iClipPyL_tup, label_tup,
                    norm_tup = None,
                    flank = 0, 
                    stranded = True,
                    cleav_rateMode = False,
                    oneD_rateMode = True,
                    rateMode_cut = 15
                ):
    
    print('Flank', flank) #debugging
    gmGraphParamsD = bed_gen.calcGraphParams(flankint = flank)
    siteTup = (  gmGraphParamsD['chrom'], 
                 gmGraphParamsD['gmIntervalStart'], 
                 gmGraphParamsD['gmIntervalEnd'], 
                 gmGraphParamsD['strand']
               )
    vectorTupL = []
    if not norm_tup: norm_tup = (None, ) * len(iClipPyL_tup)
    for iClipPyL, label, norm in zip(iClipPyL_tup, label_tup, norm_tup):
        
        #if this data is from a stranded assay then add the data
        #from the gene model strand only.
        
        #NOTE: sigTup denotes "signature tuple"
        sigTup = iClipPyL.siteTup2sigTup(siteTup, 
                                         cleav_rateMode = cleav_rateMode, 
                                         oneD_rateMode = oneD_rateMode, 
                                         stranded = stranded,
                                         rateMode_cut = rateMode_cut)
        rawCovL, cleavCovL, oneDcovL = sigTup
        
        if norm:
            rawCovL = [rawCov/norm for rawCov in rawCovL]
            if not cleav_rateMode:
                cleavCovL = [cleavCov/norm for cleavCov in cleavCovL]
        else:
            pass
        
        
        vectorTupL.append( (rawCovL, cleavCovL, oneDcovL) )
    
    
    return vectorTupL

def hitsCLIPvectors(bed_gen, 
                    hitsClipPyL_tup, label_tup,
                    norm_tup = None,
                    flank = 0, 
                    stranded = True,
                    cleav_rateMode = False,
                    oneD_rateMode = True,
                    rateMode_cut = 15
                ):
    
    print('Flank', flank) #debugging
    gmGraphParamsD = bed_gen.calcGraphParams(flankint = flank)
    siteTup = (  gmGraphParamsD['chrom'], 
                 gmGraphParamsD['gmIntervalStart'], 
                 gmGraphParamsD['gmIntervalEnd'], 
                 gmGraphParamsD['strand']
               )
    vectorTupL = []
    if not norm_tup: norm_tup = (None, ) * len(hitsClipPyL_tup)
    for hitsClipPyL, label, norm in zip(hitsClipPyL_tup, label_tup, norm_tup):
        
        #if this data is from a stranded assay then add the data
        #from the gene model strand only.
        
        #NOTE: sigTup denotes "signature tuple"
        sigTup = hitsClipPyL.siteTup2sigTup(siteTup, 
                                         cleav_rateMode = cleav_rateMode, 
                                         oneD_rateMode = oneD_rateMode, 
                                         stranded = stranded,
                                         rateMode_cut = rateMode_cut)
        rawCovL, cleavCovL, oneDcovL = sigTup
        
        if norm:
            rawCovL = [rawCov/norm for rawCov in rawCovL]
            if not cleav_rateMode:
                cleavCovL = [cleavCov/norm for cleavCov in cleavCovL]
        else:
            pass
        
        
        vectorTupL.append( (rawCovL, cleavCovL, oneDcovL) )
    
    
    return vectorTupL
