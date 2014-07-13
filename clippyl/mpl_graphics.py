import pysam
import math

from clippyl.sqlite_io import ReadidSQLite, Bed6SQLite
from clippyl.vector_factory import build_hitsclip_vectors
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.backends.backend_pdf import PdfPages

def hits_clip_plot(   bed_gen,
                      bam_fh_l,
                      readid_db_fh_l = None,
                      label_l = None,
                      norm_factor_l = None,
                      ciselement_db_fh_l = None,
                      ciselement_label_l = None,
                      ciselement_color_l = None,
                      flank = 0, 
                      all_adapter_clipped = False,
                      uniq_only = True,
                      oneD_rate_mode = True,
                      rate_cutoff = 15,
                      cleav_rate_mode = False,
                      stranded = True):
    '''
    writeme
    '''
    
    # Instantiate figure instance that will hold multiple axes objects
    fig = plt.figure(figsize=(11,8.5))
    
    # the region defined in the bed file will determined the graph limits
    # the flank argument expands the graph to include additional nts on 
    # each side
    print('Flank', flank) #debugging
    bed_d = bed_gen.calc_graph_limits(flank_size = flank)
    
    if not norm_factor_l:
        norm_factor_l = [None, ] * len(bam_fh_l)
    
    if not label_l:
        label_l = ['no label', ] * len(bam_fh_l)
    
    if not readid_db_fh_l:
        readid_db_fh_l = [None, ] * len(bam_fh_l)
    
    # these containers for max values will help determine y-axis length
    cleav_max = 5
    oneD_max = 5
    raw_cover_max = 5
    
    for bam_fh, label, norm_factor, readid_db_fh in zip(bam_fh_l,
                                                        label_l,
                                                        norm_factor_l,
                                                        readid_db_fh_l):
        
        #if this data is from a stranded assay then add the data
        #from the gene model strand only.
        
        #NOTE: sigTup denotes "signature tuple"
        # Each data chromDict will be queried for data within the geneModel bounds
        ome_coords = (bed_d['ref'],
                      bed_d['graph_start'],
                      bed_d['graph_end'],
                      bed_d['strand'])
        
        f = build_hitsclip_vectors()
        t = f( ome_coords,
               bam_fh, 
               cleaved_readid_db_conn = readid_db_fh,
               all_adapter_clipped = False,
               uniq_only = True,
               oneD_rate_mode = True,
               rate_cutoff = 15,
               cleav_rate_mode = False,
               stranded = True)
        
        # recover the coverage vectors for the region defined by ome_coords
        raw_cover_l, oneD_l, cleavage_l = t
        
        if norm_factor:
            raw_cover_l = [raw_cover/norm_factor for raw_cover in raw_cover_l]
            if not cleav_rate_mode:
                cleavage_l = [cleavage/norm_factor for cleavage in cleavage_l]
            if not oneD_rate_mode:
                oneD_l = [oneD/norm_factor for oneD in oneD_l]
        else:
            pass
        
        #### Lineplot part
        # Instantiate axes that will be populated with lineplot and legend (pyplot)
        # [left, bottom, width, height]
        #TODO: move legend to its own axes (at top)
        #TODO: remove x-axis labels from ax3 and ax2
        #TODO: widen left margin to accomodate y-axis labels
        # note: axes instantiated with pyplot in this way will be drawn from bottom 
        # such that cleav_ax is bottom and raw_ax is top
        
        cleav_ax = plt.axes([0.085, 0.15, 0.9, 0.2])
        cleav_ax.plot(cleavage_l, label = label)
        if max(cleavage_l, key = int) > cleav_max:
            cleav_max = max(cleavage_l, key = int)
        
        oneD_ax = plt.axes([0.085, 0.35, 0.9, 0.2])
        oneD_ax.plot(oneD_l, label = label)
        if max(oneD_l, key = int) > oneD_max:
            oneD_max = max(oneD_l, key = int)
        
        raw_ax = plt.axes([0.085, 0.55, 0.9, 0.2])
        raw_ax.plot(raw_cover_l, label = label)
        if max(raw_cover_l, key = int) > raw_cover_max:
            raw_cover_max = max(raw_cover_l, key = int)
    
    # delineate query region and flank region with vspan
    q_region_ax = fig.add_axes([0.085, 0.08, 0.9, 0.02], sharex = raw_ax)
    x_axis_len = bed_d['graph_end'] - bed_d['graph_start']
    q_region_ax.axvspan(flank, x_axis_len-flank, facecolor='k', alpha=0.5)
    plt.tick_params( axis='y',       # changes apply to the x-axis
                     which='both',      # both major and minor ticks are affected
                     #bottom='off',      # ticks along the bottom edge are off
                     #top='off',         # ticks along the top edge are off
                     left='off', 
                     right='off', 
                     #labelbottom='off', # labels along the bottom edge are off
                     labelleft='off',
                     labelright='off' )
    
    # configure x and y axis for the raw coverage plot
    plt.sca(raw_ax)
    plt.xlim( (0, x_axis_len) )
    plt.ylim( (0, raw_cover_max) )
    plt.ylabel( 'read coverage' )
    plt.tick_params( axis='x',          # changes apply to the x-axis
                     which='both',      # both major and minor ticks are affected
                     #bottom='off',      # ticks along the bottom edge are off
                     #top='off',         # ticks along the top edge are off
                     labelbottom='off') # labels along the bottom edge are off
    locs, labels = plt.yticks()
    plt.yticks(locs[1:])
    
    # configure x and y axis for the oneD plot
    plt.sca(oneD_ax)
    plt.xlim( (0, x_axis_len) )
    if oneD_rate_mode:
        plt.ylim( (0, 1) )
    else:
        plt.ylim( (0, oneD_max) )
    plt.ylabel( '1D rate' )
    plt.tick_params( axis='x',          # changes apply to the x-axis
                     which='both',      # both major and minor ticks are affected
                     #bottom='off',      # ticks along the bottom edge are off
                     #top='off',         # ticks along the top edge are off
                     labelbottom='off') # labels along the bottom edge are off
    locs, labels = plt.yticks()
    plt.yticks(locs[1:-1])
    
    # configure x and y axis for the cleavage plot
    plt.sca(cleav_ax)
    plt.xlim( (0, x_axis_len) )
    if cleav_rate_mode:
        plt.ylim( (0, 1) )
    else:
        plt.ylim( (0, cleav_max) )
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
    
    text_s = '\n'.join( ['Coordinates: ' + bed_d['ref'] + ':' \
                                         + str(bed_d['start']) + '-' \
                                         + str(bed_d['end']), 
                         'Name: ' + bed_d['name'], 
                         'Score: ' + str(bed_d['score']), 
                         'Strand: ' + bed_d['strand']] )
    
    plt.text( 0.7, 0.7, text_s, 
              horizontalalignment='center',
              verticalalignment='center',
              transform=legend_ax.transAxes, 
              bbox=dict() )
    
    if not ciselement_color_l:
        ciselement_color_l = ['g', ] * len(ciselement_db_fh_l)
    
    # check for motif intersect
    for db_fh, label, color in zip(ciselement_db_fh_l, 
                                   ciselement_label_l, 
                                   ciselement_color_l):
        
        l = db_fh.ome_coord_lookup( ome_coords )
        if l:
            print('cis-element intersection detected') #debugging
            for ref, start, end, name, score, strand in l:
                ce_vspan_start = start - bed_d['graph_start']
                ce_vspan_end = end - bed_d['graph_start']
                if ce_vspan_start < 0:
                    # the motif is only partially overlapping
                    ce_vspan_start = 0
                if ce_vspan_end > x_axis_len:
                    # the motif is only partially overlapping
                    ce_vspan_end = x_axis_len
                
                #make sure everything is in the same orientation
                if bed_d['strand'] == '-':
                    i = len(raw_cover_l) - ce_vspan_end
                    j = len(raw_cover_l) - ce_vspan_start
                else:
                    i = ce_vspan_start
                    j = ce_vspan_end
                raw_ax.axvspan(i, j, alpha=0.25, facecolor = color)
    
    return fig

#TODO
def par_clip_plot(   bed_gen,
                      bam_fh_l,
                      readid_db_fh_l = None,
                      label_l = None,
                      norm_factor_l = None,
                      ciselement_db_fh_l = None,
                      ciselement_label_l = None,
                      ciselement_color_l = None,
                      flank = 0, 
                      all_adapter_clipped = False,
                      uniq_only = True,
                      oneD_rate_mode = True,
                      rate_cutoff = 15,
                      cleav_rate_mode = False,
                      stranded = True):
    '''
    writeme
    '''
    
    # Instantiate figure instance that will hold multiple axes objects
    fig = plt.figure(figsize=(11,8.5))
    
    # the region defined in the bed file will determined the graph limits
    # the flank argument expands the graph to include additional nts on 
    # each side
    print('Flank', flank) #debugging
    bed_d = bed_gen.calc_graph_limits(flank_size = flank)
    
    if not norm_factor_l:
        norm_factor_l = [None, ] * len(bam_fh_l)
    
    if not label_l:
        label_l = ['no label', ] * len(bam_fh_l)
    
    if not readid_db_fh_l:
        readid_db_fh_l = [None, ] * len(bam_fh_l)
    
    # these containers for max values will help determine y-axis length
    cleav_max = 5
    oneD_max = 5
    raw_cover_max = 5
    
    for bam_fh, label, norm_factor, readid_db_fh in zip(bam_fh_l,
                                                        label_l,
                                                        norm_factor_l,
                                                        readid_db_fh_l):
        
        #if this data is from a stranded assay then add the data
        #from the gene model strand only.
        
        #NOTE: sigTup denotes "signature tuple"
        # Each data chromDict will be queried for data within the geneModel bounds
        ome_coords = (bed_d['ref'],
                      bed_d['graph_start'],
                      bed_d['graph_end'],
                      bed_d['strand'])
        
        f = build_hitsclip_vectors()
        t = f( ome_coords,
               bam_fh, 
               cleaved_readid_db_conn = readid_db_fh,
               all_adapter_clipped = False,
               uniq_only = True,
               oneD_rate_mode = True,
               rate_cutoff = 15,
               cleav_rate_mode = False,
               stranded = True)
        
        # recover the coverage vectors for the region defined by ome_coords
        raw_cover_l, oneD_l, cleavage_l = t
        
        if norm_factor:
            raw_cover_l = [raw_cover/norm_factor for raw_cover in raw_cover_l]
            if not cleav_rate_mode:
                cleavage_l = [cleavage/norm_factor for cleavage in cleavage_l]
            if not oneD_rate_mode:
                oneD_l = [oneD/norm_factor for oneD in oneD_l]
        else:
            pass
        
        #### Lineplot part
        # Instantiate axes that will be populated with lineplot and legend (pyplot)
        # [left, bottom, width, height]
        #TODO: move legend to its own axes (at top)
        #TODO: remove x-axis labels from ax3 and ax2
        #TODO: widen left margin to accomodate y-axis labels
        # note: axes instantiated with pyplot in this way will be drawn from bottom 
        # such that cleav_ax is bottom and raw_ax is top
        
        cleav_ax = plt.axes([0.085, 0.15, 0.9, 0.2])
        cleav_ax.plot(cleavage_l, label = label)
        if max(cleavage_l, key = int) > cleav_max:
            cleav_max = max(cleavage_l, key = int)
        
        oneD_ax = plt.axes([0.085, 0.35, 0.9, 0.2])
        oneD_ax.plot(oneD_l, label = label)
        if max(oneD_l, key = int) > oneD_max:
            oneD_max = max(oneD_l, key = int)
        
        raw_ax = plt.axes([0.085, 0.55, 0.9, 0.2])
        raw_ax.plot(raw_cover_l, label = label)
        if max(raw_cover_l, key = int) > raw_cover_max:
            raw_cover_max = max(raw_cover_l, key = int)
    
    # delineate query region and flank region with vspan
    q_region_ax = fig.add_axes([0.085, 0.08, 0.9, 0.02], sharex = raw_ax)
    x_axis_len = bed_d['graph_end'] - bed_d['graph_start']
    q_region_ax.axvspan(flank, x_axis_len-flank, facecolor='k', alpha=0.5)
    plt.tick_params( axis='y',       # changes apply to the x-axis
                     which='both',      # both major and minor ticks are affected
                     #bottom='off',      # ticks along the bottom edge are off
                     #top='off',         # ticks along the top edge are off
                     left='off', 
                     right='off', 
                     #labelbottom='off', # labels along the bottom edge are off
                     labelleft='off',
                     labelright='off' )
    
    # configure x and y axis for the raw coverage plot
    plt.sca(raw_ax)
    plt.xlim( (0, x_axis_len) )
    plt.ylim( (0, raw_cover_max) )
    plt.ylabel( 'read coverage' )
    plt.tick_params( axis='x',          # changes apply to the x-axis
                     which='both',      # both major and minor ticks are affected
                     #bottom='off',      # ticks along the bottom edge are off
                     #top='off',         # ticks along the top edge are off
                     labelbottom='off') # labels along the bottom edge are off
    locs, labels = plt.yticks()
    plt.yticks(locs[1:])
    
    # configure x and y axis for the oneD plot
    plt.sca(oneD_ax)
    plt.xlim( (0, x_axis_len) )
    if oneD_rate_mode:
        plt.ylim( (0, 1) )
    else:
        plt.ylim( (0, oneD_max) )
    plt.ylabel( '1D rate' )
    plt.tick_params( axis='x',          # changes apply to the x-axis
                     which='both',      # both major and minor ticks are affected
                     #bottom='off',      # ticks along the bottom edge are off
                     #top='off',         # ticks along the top edge are off
                     labelbottom='off') # labels along the bottom edge are off
    locs, labels = plt.yticks()
    plt.yticks(locs[1:-1])
    
    # configure x and y axis for the cleavage plot
    plt.sca(cleav_ax)
    plt.xlim( (0, x_axis_len) )
    if cleav_rate_mode:
        plt.ylim( (0, 1) )
    else:
        plt.ylim( (0, cleav_max) )
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
    
    text_s = '\n'.join( ['Coordinates: ' + bed_d['ref'] + ':' \
                                         + str(bed_d['start']) + '-' \
                                         + str(bed_d['end']), 
                         'Name: ' + bed_d['name'], 
                         'Score: ' + str(bed_d['score']), 
                         'Strand: ' + bed_d['strand']] )
    
    plt.text( 0.7, 0.7, text_s, 
              horizontalalignment='center',
              verticalalignment='center',
              transform=legend_ax.transAxes, 
              bbox=dict() )
    
    if not ciselement_color_l:
        ciselement_color_l = ['g', ] * len(ciselement_db_fh_l)
    
    # check for motif intersect
    for db_fh, label, color in zip(ciselement_db_fh_l, 
                                   ciselement_label_l, 
                                   ciselement_color_l):
        
        l = db_fh.ome_coord_lookup( ome_coords )
        if l:
            print('cis-element intersection detected') #debugging
            for ref, start, end, name, score, strand in l:
                ce_vspan_start = start - bed_d['graph_start']
                ce_vspan_end = end - bed_d['graph_start']
                if ce_vspan_start < 0:
                    # the motif is only partially overlapping
                    ce_vspan_start = 0
                if ce_vspan_end > x_axis_len:
                    # the motif is only partially overlapping
                    ce_vspan_end = x_axis_len
                
                #make sure everything is in the same orientation
                if bed_d['strand'] == '-':
                    i = len(raw_cover_l) - ce_vspan_end
                    j = len(raw_cover_l) - ce_vspan_start
                else:
                    i = ce_vspan_start
                    j = ce_vspan_end
                raw_ax.axvspan(i, j, alpha=0.25, facecolor = color)
    
    return fig


#################################LEGACY CODE FROM HERE DOWN KEPT FOR POSTERITY (will delete after pysam migration is complete)

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
        # note: axes instantiated with pyplot in this way will be drawn from bottom 
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

