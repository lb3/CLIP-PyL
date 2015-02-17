import os

from clippyl.build_readid_db import build_ReadidSQLite_dbs
from clippyl.coverage_graphics import hitsclip_graphics
from clippyl.sample_data.paths import hitsclip_discardUnclipped_fq_dir

if __name__ == '__main__':
    # getting relevant sample data filepaths
    fp_l = sorted(os.listdir(hitsclip_discardUnclipped_fq_dir()))
    fp_l = [os.path.join(hitsclip_discardUnclipped_fq_dir(), fp) for fp in fp_l]
    fq_fp_l = [fp for fp in fp_l if fp.split('.')[-2:] == ['fq','gz']]
    #print (fq_fp_l) #debugging
    #NOTE: expected sort order is:
    #SLBP_CLIP_high_MW_high_MNase.HISTONLY.discardUnclipped.fq.gz
    #SLBP_CLIP_high_MW_low_MNase.HISTONLY.discardUnclipped.fq.gz
    #SLBP_CLIP_low_MW_high_MNase.HISTONLY.discardUnclipped.fq.gz
    #SLBP_CLIP_low_MW_low_MNase.HISTONLY.discardUnclipped.fq.gz
    #mock_CLIP_high_MW_high_MNase.HISTONLY.discardUnclipped.fq.gz
    #mock_CLIP_low_MW_high_MNase.HISTONLY.discardUnclipped.fq.gz
    
    build_ReadidSQLite_dbs(fq_fp_l)
    
    #TODO: parameterize the normalization factor because histonly files do not contain all reads
    #number of mapped reads from full bam file are supplied in this 
    #parameter because HISTONLY files do not contain all reads:
    n_mapped_reads_l = [17931502, #SLBP_CLIP_high_MW_high_MNase.bam
                        14698879, #SLBP_CLIP_high_MW_low_MNase.bam
                        20020463, #SLBP_CLIP_low_MW_high_MNase.bam
                        15256431, #SLBP_CLIP_low_MW_low_MNase.bam
                        18992515, #mock_CLIP_high_MW_high_MNase.bam
                        21083226, #mock_CLIP_low_MW_high_MNase.bam
                       ]
    
    hitsclip_graphics(
