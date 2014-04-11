

#TODO: add support for bed12 for spliced rna that
# uses the siteTup2sigTup methods and then combines
# the splices the coverageL together, marks the union 
# in the plot
# example: /storage/Ziggy_BigGuy/LB_Bioinformatics/data_Projects/Brooks_HITS_CLIP_SLBP/chromPlots/LB_draw_hist_gene_model_linegraph_vRawCover.py

import math

import matplotlib.pyplot as plt
from matplotlib import ticker

from clipPyL.flatfile_parsing import Bed6Reader

#this function takes a tup of clipPyL db that hav been loaded
#NOTE THIS CAN ALSO BE USED TO PLOT iCLIP data (same signature components)
def hitsCliPyL_plot(   bed_gen, 
                        hitsClipPyL_tup, label_tup,
                        motif_db_tup = (), motif_label_tup = (),
                        norm_tup = None,
                        motif_color_tup = None,
                        flank = 0, 
                        yLabel = 'Y-Axis Label', 
                        graphTitle = 'Title', 
                        yAx = None, 
                        stranded = True,
                        oneD_rateMode = True,
                        cleav_rateMode = False,
                        rateMode_cut = 15
                       ):
    '''
    The plot is calculated for the region of the chromosome specified by the gene model
    graph parameters stored in gmGraphParamsD.
    NOTE: currently the only method that produces such a dict is calcGMgraphParams method
    of the Bed12Reader and there is the same method in Bed6Reader.
    
    The data encoded in the "data chromDict" files specified in the tuple tupOfDataChromDFP
    will be plotted as linegraphs and the legend with show the labels in tupOfDataLabels.
    
    Additionally, the motif locations spanned by the "motif chromDict" files declared in 
    the tuple tupOfMotifChromDFP will be marked in the line graphs; if any are present within
    the gene model region then vspan will be used to mark them.
    
    Optionally:
    You may declare the y axis label and the graph title.
    You may declare the color tuple used with multiple motif chromDicts. If the length
    of the color tuple is shorter than the tuple of motif chromDicts then it is recycled.
    You may declare the yAx size with the yAx arg
    
    NOTE: all gene models are plotted with a left-right orientation; negative strand data is reversed
    '''
    
    #NOTE: template code at /storage/Ziggy_BigGuy/LB_Bioinformatics/LB_Scripts/BioPyL/BioPyL-0.1dev/incorporateme/LB_Plotter_Master_v16/LB_Coverage_Plotterator_Legacy_Code/Spyder_Project_Coverage_Plotterator_v00/Plotterator_Master_v00.py
    
    #TODO: allow tuple of colors to be specified for the data chromDicts line plot lines
    #TODO: allow tuple of colors to be specified for the motif chromDicts span regions
    #TODO: make legend for motif chromDict labels
    #TODO: make an "exons only" argument option
    #TODO: make positive AND negative strand plot windows as in latest iteration...may want to build that as and entirely different plotting function
    
    print('Flank', flank) #debugging
    gmGraphParamsD = bed_gen.calcGraphParams(flankint = flank)
    siteTup = (  gmGraphParamsD['chrom'], 
                 gmGraphParamsD['gmIntervalStart'], 
                 gmGraphParamsD['gmIntervalEnd'], 
                 gmGraphParamsD['strand']
               )
    
    
    # Instantiate figure instance that will hold multiple axes objects
    fig = plt.figure(figsize=(11,8.5))
    
    xAxMax = (gmGraphParamsD['gmIntervalEnd']) - (gmGraphParamsD['gmIntervalStart'])
    print('xAxMax',xAxMax) #debugging
    
    yAxMax = 10
    
    rawMax = 0
    cleavMax = 1
    oneDmax = 1
    
    if not norm_tup: norm_tup = (None, ) * len(hitsClipPyL_tup)
    for hitsClipPyL, label, norm in zip(hitsClipPyL_tup, label_tup, norm_tup):
        
        #if this data is from a stranded assay then add the data
        #from the gene model strand only.
        
        #NOTE: sigTup denotes "signature tuple"
        # Each data chromDict will be queried for data within the geneModel bounds
        sigTup = hitsClipPyL.siteTup2sigTup(siteTup, 
                                            oneD_rateMode = oneD_rateMode, 
                                            cleav_rateMode = cleav_rateMode, 
                                            stranded = stranded,
                                            rateMode_cut = rateMode_cut)
        rawCovL, cleavCovL, oneDcovL = sigTup
        
        if norm:
            rawCovL = [rawCov/norm for rawCov in rawCovL]
            if not cleav_rateMode:
                cleavCovL = [cleavCov/norm for cleavCov in cleavCovL]
        else:
            pass
        
        #### Lineplot part
        # Instantiate axes that will be populated with lineplot and legend (pyplot)
        # [left, bottom, width, height]
        #TODO: move legend to its own axes (at top)
        #TODO: remove x-axis labels from ax3 and ax2
        #TODO: widen left margin to accomodate y-axis labels
        # note: axes instatiated with pyplot in this way will be drawn from bottom 
        # such that cleav_ax is bottom and raw_ax is top
        
        cleav_ax = plt.axes([0.085, 0.15, 0.9, 0.2])
        cleav_ax.plot(cleavCovL, label = label)
        if max(cleavCovL) > cleavMax: cleavMax = max(cleavCovL)
        
        oneD_ax = plt.axes([0.085, 0.35, 0.9, 0.2])
        oneD_ax.plot(oneDcovL, label = label)
        if max(oneDcovL) > oneDmax: oneDmax = max(oneDcovL)
        
        raw_ax = plt.axes([0.085, 0.55, 0.9, 0.2])
        raw_ax.plot(rawCovL, label = label)
        #store max data value to use for explicit declaration of yAxis height downstream
        if max(rawCovL) > rawMax: rawMax = max(rawCovL)
        
        if yAxMax < max(rawCovL, key=int): yAxMax = max(rawCovL, key=int) + 1
    
    ##mark genemodel region with vspan
    #TODO: use matplotlib.pyplot.broken_barh or revert to colorbar! for gene model ranges
    gene_model_ax = fig.add_axes([0.085, 0.08, 0.9, 0.02], sharex = raw_ax)
    gene_model_ax.axvspan(flank, xAxMax-flank, facecolor='k', alpha=0.5)
    plt.tick_params( axis='y',       # changes apply to the x-axis
                     which='both',      # both major and minor ticks are affected
                     #bottom='off',      # ticks along the bottom edge are off
                     #top='off',         # ticks along the top edge are off
                     left='off', 
                     right='off', 
                     #labelbottom='off', # labels along the bottom edge are off
                     labelleft='off',
                     labelright='off' )
    
    # create a list to explicitly set the min and max of the x and y axes, 
    # such that v = [xmin, xmax, ymin, ymax].
    # for more options see pyplot api documentation; 
    # http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.axis
    if yAx != None: 
        v = [0, xAxMax, 0, yAx]
    else:
        v = [0, xAxMax, 0, yAxMax]
    
    plt.sca(raw_ax)
    plt.axis(v) #note: this call to axis() must follow the plot() call
    plt.ylabel( 'read coverage' )
    plt.tick_params( axis='x',          # changes apply to the x-axis
                     which='both',      # both major and minor ticks are affected
                     #bottom='off',      # ticks along the bottom edge are off
                     #top='off',         # ticks along the top edge are off
                     labelbottom='off') # labels along the bottom edge are off
    locs, labels = plt.yticks()
    plt.yticks(locs[1:])
    
    plt.sca(oneD_ax)
    plt.xlim( (0, xAxMax) )
    plt.ylim( (0, oneDmax) )
    plt.ylabel( '1D rate' )
    plt.tick_params( axis='x',          # changes apply to the x-axis
                     which='both',      # both major and minor ticks are affected
                     #bottom='off',      # ticks along the bottom edge are off
                     #top='off',         # ticks along the top edge are off
                     labelbottom='off') # labels along the bottom edge are off
    locs, labels = plt.yticks()
    plt.yticks(locs[1:-1])
    
    plt.sca(cleav_ax)
    plt.xlim( (0, xAxMax) )
    plt.ylim( (0, cleavMax) )
    plt.ylabel( 'fragment termini' )
    locs, labels = plt.yticks()
    plt.yticks(locs[1:-1]) #remove max and min ticks so that they don't overlap with adjacent plots on the display
    
    # this will be an axis area to draw the legend and some
    # text info about the gene model
    legend_ax = plt.axes([0.085, 0.75, 0.9, 0.2], frameon=False)
    handles, labels = raw_ax.get_legend_handles_labels()
    ncol = math.ceil( len(labels) / 4 )
    legend_ax.legend(handles, labels, ncol = ncol, loc=2, fontsize=12)
    plt.tick_params( axis='both',       # changes apply to the x-axis
                     which='both',      # both major and minor ticks are affected
                     bottom='off',      # ticks along the bottom edge are off
                     top='off',         # ticks along the top edge are off
                     left='off', 
                     right='off', 
                     labelbottom='off', # labels along the bottom edge are off
                     labelleft='off',
                     labelright='off' )
    
    text_s = '\n'.join( ['Coordinates: ' + bed_gen.d['chrom'] + ':' \
                                         + str(bed_gen.d['chromStart']) + '-' \
                                         + str(bed_gen.d['chromEnd']), 
                         'Name: ' + bed_gen.d['name'], 
                         'Score: ' + str(bed_gen.d['score']), 
                         'Strand: ' + bed_gen.d['strand']] )
    
    plt.text( 0.7, 0.7, text_s, 
              horizontalalignment='center',
              verticalalignment='center',
              transform=legend_ax.transAxes, 
              bbox=dict() )
    
    if not motif_color_tup: motif_color_tup = ('g', ) * len(motif_db_tup)
    # check for motif intersect
    for motif_db, motif_label, motif_color in zip(motif_db_tup, motif_label_tup, motif_color_tup):
        
        print(siteTup)
        motif_siteTup_list = motif_db.siteTup_lookup(siteTup)
        if motif_siteTup_list:
            print('motif intersection detected')
            for chrom, chromStart, chromEnd, name, score, strand in motif_siteTup_list:
                motif_start = chromStart - gmGraphParamsD['gmIntervalStart']
                motif_end = chromEnd - gmGraphParamsD['gmIntervalStart']
                if motif_start < 0:
                    # if the motif is only partially overlapping
                    motif_start = 0
                #make sure everything is in the same orientation
                if gmGraphParamsD['strand'] == '-':
                    i = len(rawCovL) - motif_end
                    j = len(rawCovL) - motif_start
                else:
                    i = motif_start
                    j = motif_end
                raw_ax.axvspan(i, j, alpha=0.25, facecolor = motif_color)
    
    return fig

#this function takes a tup of clipPyL db that hav been loaded
def iCliPyL_plot(   bed_gen, 
                        iClipPyL_tup, label_tup,
                        motif_db_tup = (), motif_label_tup = (),
                        norm_tup = None,
                        motif_color_tup = None,
                        flank = 0, 
                        yLabel = 'Y-Axis Label', 
                        graphTitle = 'Title', 
                        yAx = None, 
                        stranded = True,
                        oneD_rateMode = True,
                        cleav_rateMode = False
                       ):
    '''
    The plot is calculated for the region of the chromosome specified by the gene model
    graph parameters stored in gmGraphParamsD.
    NOTE: currently the only method that produces such a dict is calcGMgraphParams method
    of the Bed12Reader and there is the same method in Bed6Reader.
    
    The data encoded in the "data chromDict" files specified in the tuple tupOfDataChromDFP
    will be plotted as linegraphs and the legend with show the labels in tupOfDataLabels.
    
    Additionally, the motif locations spanned by the "motif chromDict" files declared in 
    the tuple tupOfMotifChromDFP will be marked in the line graphs; if any are present within
    the gene model region then vspan will be used to mark them.
    
    Optionally:
    You may declare the y axis label and the graph title.
    You may declare the color tuple used with multiple motif chromDicts. If the length
    of the color tuple is shorter than the tuple of motif chromDicts then it is recycled.
    You may declare the yAx size with the yAx arg
    
    NOTE: all gene models are plotted with a left-right orientation; negative strand data is reversed
    '''
    
    #NOTE: template code at /storage/Ziggy_BigGuy/LB_Bioinformatics/LB_Scripts/BioPyL/BioPyL-0.1dev/incorporateme/LB_Plotter_Master_v16/LB_Coverage_Plotterator_Legacy_Code/Spyder_Project_Coverage_Plotterator_v00/Plotterator_Master_v00.py
    
    #TODO: allow tuple of colors to be specified for the data chromDicts line plot lines
    #TODO: allow tuple of colors to be specified for the motif chromDicts span regions
    #TODO: make legend for motif chromDict labels
    #TODO: make an "exons only" argument option
    #TODO: make positive AND negative strand plot windows as in latest iteration...may want to build that as and entirely different plotting function
    
    print('Flank', flank) #debugging
    gmGraphParamsD = bed_gen.calcGraphParams(flankint = flank)
    siteTup = (  gmGraphParamsD['chrom'], 
                 gmGraphParamsD['gmIntervalStart'], 
                 gmGraphParamsD['gmIntervalEnd'], 
                 gmGraphParamsD['strand']
               )
    
    
    # Instantiate figure instance that will hold multiple axes objects
    fig = plt.figure(figsize=(11,8.5))
    
    xAxMax = (gmGraphParamsD['gmIntervalEnd']) - (gmGraphParamsD['gmIntervalStart'])
    print('xAxMax',xAxMax) #debugging
    
    yAxMax = 10
    
    rawMax = 0
    cleavMax = 1
    oneDmax = 1
    
    if not norm_tup: norm_tup = (None, ) * len(iClipPyL_tup)
    for iClipPyL, label, norm in zip(iClipPyL_tup, label_tup, norm_tup):
        
        #if this data is from a stranded assay then add the data
        #from the gene model strand only.
        
        #NOTE: sigTup denotes "signature tuple"
        # Each data chromDict will be queried for data within the geneModel bounds
        sigTup = iClipPyL.siteTup2sigTup(siteTup, 
                                         oneD_rateMode = oneD_rateMode, 
                                         cleav_rateMode = cleav_rateMode, 
                                         stranded = stranded)
        rawCovL, cleavCovL, oneDcovL = sigTup
        
        if norm:
            rawCovL = [rawCov/norm for rawCov in rawCovL]
            if not cleav_rateMode:
                cleavCovL = [cleavCov/norm for cleavCov in cleavCovL]
        else:
            pass
        
        #### Lineplot part
        # Instantiate axes that will be populated with lineplot and legend (pyplot)
        # [left, bottom, width, height]
        #TODO: move legend to its own axes (at top)
        #TODO: remove x-axis labels from ax3 and ax2
        #TODO: widen left margin to accomodate y-axis labels
        # note: axes instatiated with pyplot in this way will be drawn from bottom 
        # such that cleav_ax is bottom and raw_ax is top
        
        cleav_ax = plt.axes([0.085, 0.15, 0.9, 0.2])
        cleav_ax.plot(cleavCovL, label = label)
        if max(cleavCovL) > cleavMax: cleavMax = max(cleavCovL)
        
        oneD_ax = plt.axes([0.085, 0.35, 0.9, 0.2])
        oneD_ax.plot(oneDcovL, label = label)
        if max(oneDcovL) > oneDmax: oneDmax = max(oneDcovL)
        
        raw_ax = plt.axes([0.085, 0.55, 0.9, 0.2])
        raw_ax.plot(rawCovL, label = label)
        #store max data value to use for explicit declaration of yAxis height downstream
        if max(rawCovL) > rawMax: rawMax = max(rawCovL)
        
        if yAxMax < max(rawCovL, key=int): yAxMax = max(rawCovL, key=int) + 1
    
    ##mark genemodel region with vspan
    #TODO: use matplotlib.pyplot.broken_barh or revert to colorbar! for gene model ranges
    gene_model_ax = fig.add_axes([0.085, 0.08, 0.9, 0.02], sharex = raw_ax)
    gene_model_ax.axvspan(flank, xAxMax-flank, facecolor='k', alpha=0.5)
    plt.tick_params( axis='y',       # changes apply to the x-axis
                     which='both',      # both major and minor ticks are affected
                     #bottom='off',      # ticks along the bottom edge are off
                     #top='off',         # ticks along the top edge are off
                     left='off', 
                     right='off', 
                     #labelbottom='off', # labels along the bottom edge are off
                     labelleft='off',
                     labelright='off' )
    
    # create a list to explicitly set the min and max of the x and y axes, 
    # such that v = [xmin, xmax, ymin, ymax].
    # for more options see pyplot api documentation; 
    # http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.axis
    if yAx != None: 
        v = [0, xAxMax, 0, yAx]
    else:
        v = [0, xAxMax, 0, yAxMax]
    
    plt.sca(raw_ax)
    plt.axis(v) #note: this call to axis() must follow the plot() call
    plt.ylabel( 'read coverage' )
    plt.tick_params( axis='x',          # changes apply to the x-axis
                     which='both',      # both major and minor ticks are affected
                     #bottom='off',      # ticks along the bottom edge are off
                     #top='off',         # ticks along the top edge are off
                     labelbottom='off') # labels along the bottom edge are off
    locs, labels = plt.yticks()
    plt.yticks(locs[1:])
    
    plt.sca(oneD_ax)
    plt.xlim( (0, xAxMax) )
    plt.ylim( (0, oneDmax) )
    plt.ylabel( '1D rate' )
    plt.tick_params( axis='x',          # changes apply to the x-axis
                     which='both',      # both major and minor ticks are affected
                     #bottom='off',      # ticks along the bottom edge are off
                     #top='off',         # ticks along the top edge are off
                     labelbottom='off') # labels along the bottom edge are off
    locs, labels = plt.yticks()
    plt.yticks(locs[1:-1])
    
    plt.sca(cleav_ax)
    plt.xlim( (0, xAxMax) )
    plt.ylim( (0, cleavMax) )
    plt.ylabel( 'RT termini' )
    locs, labels = plt.yticks()
    plt.yticks(locs[1:-1]) #remove max and min ticks so that they don't overlap with adjacent plots on the display
    
    # this will be an axis area to draw the legend and some
    # text info about the gene model
    legend_ax = plt.axes([0.085, 0.75, 0.9, 0.2], frameon=False)
    handles, labels = raw_ax.get_legend_handles_labels()
    ncol = math.ceil( len(labels) / 4 )
    legend_ax.legend(handles, labels, ncol = ncol, loc=2, fontsize=12)
    plt.tick_params( axis='both',       # changes apply to the x-axis
                     which='both',      # both major and minor ticks are affected
                     bottom='off',      # ticks along the bottom edge are off
                     top='off',         # ticks along the top edge are off
                     left='off', 
                     right='off', 
                     labelbottom='off', # labels along the bottom edge are off
                     labelleft='off',
                     labelright='off' )
    
    text_s = '\n'.join( ['Coordinates: ' + bed_gen.d['chrom'] + ':' \
                                         + str(bed_gen.d['chromStart']) + '-' \
                                         + str(bed_gen.d['chromEnd']), 
                         'Name: ' + bed_gen.d['name'], 
                         'Score: ' + str(bed_gen.d['score']), 
                         'Strand: ' + bed_gen.d['strand']] )
    
    plt.text( 0.7, 0.7, text_s, 
              horizontalalignment='center',
              verticalalignment='center',
              transform=legend_ax.transAxes, 
              bbox=dict() )
    
    if not motif_color_tup: motif_color_tup = ('g', ) * len(motif_db_tup)
    # check for motif intersect
    for motif_db, motif_label, motif_color in zip(motif_db_tup, motif_label_tup, motif_color_tup):
        
        print(siteTup)
        motif_siteTup_list = motif_db.siteTup_lookup(siteTup)
        if motif_siteTup_list:
            print('motif intersection detected')
            for chrom, chromStart, chromEnd, name, score, strand in motif_siteTup_list:
                motif_start = chromStart - gmGraphParamsD['gmIntervalStart']
                motif_end = chromEnd - gmGraphParamsD['gmIntervalStart']
                if motif_start < 0:
                    # if the motif is only partially overlapping
                    motif_start = 0
                #make sure everything is in the same orientation
                if gmGraphParamsD['strand'] == '-':
                    i = len(rawCovL) - motif_end
                    j = len(rawCovL) - motif_start
                else:
                    i = motif_start
                    j = motif_end
                raw_ax.axvspan(i, j, alpha=0.25, facecolor = motif_color)
    
    return fig

#this function takes a tup of clipPyL db that hav been loaded
def parCliPyL_plot(   bed_gen, 
                        parClipPyL_tup, label_tup,
                        motif_db_tup = (), motif_label_tup = (),
                        norm_tup = None,
                        motif_color_tup = None,
                        flank = 0, 
                        yLabel = 'Y-Axis Label', 
                        graphTitle = 'Title', 
                        yAx = None, 
                        stranded = True,
                        t_to_c_rateMode = True,
                        cleav_rateMode = False
                       ):
    '''
    The plot is calculated for the region of the chromosome specified by the gene model
    graph parameters stored in gmGraphParamsD.
    NOTE: currently the only method that produces such a dict is calcGMgraphParams method
    of the Bed12Reader and there is the same method in Bed6Reader.
    
    The data encoded in the "data chromDict" files specified in the tuple tupOfDataChromDFP
    will be plotted as linegraphs and the legend with show the labels in tupOfDataLabels.
    
    Additionally, the motif locations spanned by the "motif chromDict" files declared in 
    the tuple tupOfMotifChromDFP will be marked in the line graphs; if any are present within
    the gene model region then vspan will be used to mark them.
    
    Optionally:
    You may declare the y axis label and the graph title.
    You may declare the color tuple used with multiple motif chromDicts. If the length
    of the color tuple is shorter than the tuple of motif chromDicts then it is recycled.
    You may declare the yAx size with the yAx arg
    
    NOTE: all gene models are plotted with a left-right orientation; negative strand data is reversed
    '''
    
    #NOTE: template code at /storage/Ziggy_BigGuy/LB_Bioinformatics/LB_Scripts/BioPyL/BioPyL-0.1dev/incorporateme/LB_Plotter_Master_v16/LB_Coverage_Plotterator_Legacy_Code/Spyder_Project_Coverage_Plotterator_v00/Plotterator_Master_v00.py
    
    #TODO: allow tuple of colors to be specified for the data chromDicts line plot lines
    #TODO: allow tuple of colors to be specified for the motif chromDicts span regions
    #TODO: make legend for motif chromDict labels
    #TODO: make an "exons only" argument option
    #TODO: make positive AND negative strand plot windows as in latest iteration...may want to build that as and entirely different plotting function
    
    print('Flank', flank) #debugging
    gmGraphParamsD = bed_gen.calcGraphParams(flankint = flank)
    
    # Instantiate figure instance that will hold multiple axes objects
    fig = plt.figure(figsize=(11,8.5))
    
    xAxMax = (gmGraphParamsD['gmIntervalEnd']) - (gmGraphParamsD['gmIntervalStart'])
    print('xAxMax',xAxMax) #debugging
    
    yAxMax = 10
    
    rawMax = 0
    cleavMax = 1
    t_to_c_max = 1
    
    if not norm_tup: norm_tup = (None, ) * len(parClipPyL_tup)
    for parClipPyL, label, norm in zip(parClipPyL_tup, label_tup, norm_tup):
        
        #if this data is from a stranded assay then add the data
        #from the gene model strand only.
        
        #NOTE: sigTup denotes "signature tuple"
        # Each data chromDict will be queried for data within the geneModel bounds
        siteTup = (  gmGraphParamsD['chrom'], 
                     gmGraphParamsD['gmIntervalStart'], 
                     gmGraphParamsD['gmIntervalEnd'], 
                     gmGraphParamsD['strand']
                   )
        sigTup = parClipPyL.siteTup2sigTup(siteTup, 
                                           cleav_rateMode = cleav_rateMode, 
                                           t_to_c_rateMode = t_to_c_rateMode, 
                                           stranded = stranded)
        rawCovL, cleavCovL, t_to_c_covL = sigTup
        
        if norm:
            rawCovL = [rawCov/norm for rawCov in rawCovL]
            if not cleav_rateMode:
                cleavCovL = [cleavCov/norm for cleavCov in cleavCovL]
        else:
            pass
        
        #### Lineplot part
        # Instantiate axes that will be populated with lineplot and legend (pyplot)
        # [left, bottom, width, height]
        #TODO: move legend to its own axes (at top)
        #TODO: remove x-axis labels from ax3 and ax2
        #TODO: widen left margin to accomodate y-axis labels
        # note: axes instatiated with pyplot in this way will be drawn from bottom 
        # such that cleav_ax is bottom and raw_ax is top
        
        cleav_ax = plt.axes([0.085, 0.15, 0.9, 0.2])
        cleav_ax.plot(cleavCovL, label = label)
        if max(cleavCovL) > cleavMax: cleavMax = max(cleavCovL)
        
        t_to_c_ax = plt.axes([0.085, 0.35, 0.9, 0.2])
        t_to_c_ax.plot(t_to_c_covL, label = label)
        if max(t_to_c_covL) > t_to_c_max: t_to_c_max = max(t_to_c_covL)
        
        raw_ax = plt.axes([0.085, 0.55, 0.9, 0.2])
        raw_ax.plot(rawCovL, label = label)
        #store max data value to use for explicit declaration of yAxis height downstream
        if max(rawCovL) > rawMax: rawMax = max(rawCovL)
        
        if yAxMax < max(rawCovL, key=int): yAxMax = max(rawCovL, key=int) + 1
    
    ##mark genemodel region with vspan
    #TODO: use matplotlib.pyplot.broken_barh or revert to colorbar! for gene model ranges
    gene_model_ax = fig.add_axes([0.085, 0.08, 0.9, 0.02], sharex = raw_ax)
    gene_model_ax.axvspan(flank, xAxMax-flank, facecolor='k', alpha=0.5)
    plt.tick_params( axis='y',       # changes apply to the x-axis
                     which='both',      # both major and minor ticks are affected
                     #bottom='off',      # ticks along the bottom edge are off
                     #top='off',         # ticks along the top edge are off
                     left='off', 
                     right='off', 
                     #labelbottom='off', # labels along the bottom edge are off
                     labelleft='off',
                     labelright='off' )
    
    # create a list to explicitly set the min and max of the x and y axes, 
    # such that v = [xmin, xmax, ymin, ymax].
    # for more options see pyplot api documentation; 
    # http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.axis
    if yAx != None: 
        v = [0, xAxMax, 0, yAx]
    else:
        v = [0, xAxMax, 0, yAxMax]
    
    plt.sca(raw_ax)
    plt.axis(v) #note: this call to axis() must follow the plot() call
    plt.ylabel( 'read coverage' )
    plt.tick_params( axis='x',          # changes apply to the x-axis
                     which='both',      # both major and minor ticks are affected
                     #bottom='off',      # ticks along the bottom edge are off
                     #top='off',         # ticks along the top edge are off
                     labelbottom='off') # labels along the bottom edge are off
    locs, labels = plt.yticks()
    plt.yticks(locs[1:])
    
    plt.sca(t_to_c_ax)
    plt.xlim( (0, xAxMax) )
    plt.ylim( (0, t_to_c_max) )
    plt.ylabel( 'T->C rate' )
    plt.tick_params( axis='x',          # changes apply to the x-axis
                     which='both',      # both major and minor ticks are affected
                     #bottom='off',      # ticks along the bottom edge are off
                     #top='off',         # ticks along the top edge are off
                     labelbottom='off') # labels along the bottom edge are off
    locs, labels = plt.yticks()
    plt.yticks(locs[1:-1])
    
    plt.sca(cleav_ax)
    plt.xlim( (0, xAxMax) )
    plt.ylim( (0, cleavMax) )
    plt.ylabel( 'fragment termini' )
    locs, labels = plt.yticks()
    plt.yticks(locs[1:-1]) #remove max and min ticks so that they don't overlap with adjacent plots on the display
    
    # this will be an axis area to draw the legend and some
    # text info about the gene model
    legend_ax = plt.axes([0.085, 0.75, 0.9, 0.2], frameon=False)
    handles, labels = raw_ax.get_legend_handles_labels()
    ncol = math.ceil( len(labels) / 4 )
    legend_ax.legend(handles, labels, ncol = ncol, loc=2, fontsize=12)
    plt.tick_params( axis='both',       # changes apply to the x-axis
                     which='both',      # both major and minor ticks are affected
                     bottom='off',      # ticks along the bottom edge are off
                     top='off',         # ticks along the top edge are off
                     left='off', 
                     right='off', 
                     labelbottom='off', # labels along the bottom edge are off
                     labelleft='off',
                     labelright='off' )
    
    text_s = '\n'.join( ['Coordinates: ' + bed_gen.d['chrom'] + ':' \
                                         + str(bed_gen.d['chromStart']) + '-' \
                                         + str(bed_gen.d['chromEnd']), 
                         'Name: ' + bed_gen.d['name'], 
                         'Score: ' + str(bed_gen.d['score']), 
                         'Strand: ' + bed_gen.d['strand']] )
    
    plt.text( 0.7, 0.7, text_s, 
              horizontalalignment='center',
              verticalalignment='center',
              transform=legend_ax.transAxes, 
              bbox=dict() )
    
    if not motif_color_tup: motif_color_tup = ('g', ) * len(motif_db_tup)
    # check for motif intersect
    for motif_db, motif_label, motif_color in zip(motif_db_tup, motif_label_tup, motif_color_tup):
        
        print(siteTup)
        motif_siteTup_list = motif_db.siteTup_lookup(siteTup)
        if motif_siteTup_list:
            print('motif intersection detected')
            for chrom, chromStart, chromEnd, name, score, strand in motif_siteTup_list:
                motif_start = chromStart - gmGraphParamsD['gmIntervalStart']
                motif_end = chromEnd - gmGraphParamsD['gmIntervalStart']
                if motif_start < 0:
                    # if the motif is only partially overlapping
                    motif_start = 0
                #make sure everything is in the same orientation
                if gmGraphParamsD['strand'] == '-':
                    i = len(rawCovL) - motif_end
                    j = len(rawCovL) - motif_start
                else:
                    i = motif_start
                    j = motif_end
                raw_ax.axvspan(i, j, alpha=0.25, facecolor = motif_color)
    
    return fig

