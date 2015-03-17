import pysam
import math

from matplotlib import rcParams
import matplotlib.pyplot as plt
#import matplotlib.font_manager as mpl_font
from matplotlib import ticker
from matplotlib.backends.backend_pdf import PdfPages

from clippyl.sqlite_io import ReadidSQLite, Bed6SQLite
from clippyl.vector_factory import build_hitsclip_vectors

rcParams['font.size'] = 12
rcParams['legend.fontsize'] = 12
#alternative font configuration option:
#http://matplotlib.org/api/font_manager_api.html

def hits_clip_plot(   bed_gen,
                      bam_fh_l,
                      readid_db_fh_l = None,
                      label_l = None,
                      norm_factor_l = None,
                      ciselement_db_fh_l = None,
                      ciselement_label_l = None,
                      ciselement_color_l = None,
                      flank = 0, 
                      uniq_only = True,
                      oneD_rate_mode = True,
                      rate_cutoff = 15,
                      cleav_rate_mode = False,
                      stranded = True):
    
    # Instantiate figure instance that will hold multiple axes objects
    fig = plt.figure(figsize=(11,8.5))
    # axes objects that will be instantiated:
    #cleav_ax
    #oneD_ax
    #raw_ax
    #q_region_ax
    #legend_ax
    
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
        # note: axes instantiated with pyplot in this way will be drawn from bottom 
        # such that cleav_ax is bottom and raw_ax is top
        l = 0.1
        b = 0.15
        w = 0.85
        h = 0.2
        cleav_ax = plt.axes([l, 0.15, w, h])
        cleav_ax.plot(cleavage_l, label = label)
        if max(cleavage_l, key = int) > cleav_max:
            cleav_max = max(cleavage_l, key = int)
        
        oneD_ax = plt.axes([l, 0.35, w, h])
        oneD_ax.plot(oneD_l, label = label)
        if max(oneD_l, key = int) > oneD_max:
            oneD_max = max(oneD_l, key = int)
        
        raw_ax = plt.axes([l, 0.55, w, h])
        raw_ax.plot(raw_cover_l, label = label)
        if max(raw_cover_l, key = int) > raw_cover_max:
            raw_cover_max = max(raw_cover_l, key = int)
    
    # delineate query region and flank region with vspan
    q_region_ax = fig.add_axes([l, 0.08, w, 0.02], sharex = raw_ax)
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
    legend_ax = plt.axes([0.085, 0.75, w, h], frameon=False)
    handles, labels = raw_ax.get_legend_handles_labels()
    ncol = math.ceil( len(labels) / 4 )
    legend_ax.legend(handles, labels, ncol = ncol, loc=4)
    plt.tick_params( axis='both',       # changes apply to the x-axis
                     which='both',      # both major and minor ticks are affected
                     bottom='off',      # ticks along the bottom edge are off
                     top='off',         # ticks along the top edge are off
                     left='off', 
                     right='off', 
                     labelbottom='off', # labels along the bottom edge are off
                     labelleft='off',
                     labelright='off' )
    
    text_s = '; '.join( ['Name: ' + bed_d['name'], 
                         'Coordinates: ' + bed_d['ref'] + ':' \
                                         + str(bed_d['start']) + '-' \
                                         + str(bed_d['end']), 
                         'Strand: ' + bed_d['strand'],
                         'Score: ' + str(bed_d['score'])] )
    
    plt.text( 0, 1, text_s, 
              horizontalalignment='left',
              verticalalignment='center',
              transform=legend_ax.transAxes, 
              bbox=dict() )
    
    if ciselement_db_fh_l:
        
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

