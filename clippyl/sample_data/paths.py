"""This module contains numerous function that produce 
relative paths to the sample data"""
import os

def get_data_dir():
    data_dir = os.path.dirname(os.path.abspath(__file__))
    return data_dir

def histone_gene_bed_fp():
    fp = os.path.join(get_data_dir(),
                      'genomic_interval_queries', 
                      'RD_Histone_Genes.bed')
    return fp

def histone_sl_bed_fp():
    fp = os.path.join(get_data_dir(), 
                      'cis_element_annotations', 
                      'histone_stemloop_elements.bed')
    return fp

def histone_sl_bed_sl3_fp():
    fp = os.path.join(get_data_dir(), 
                      'cis_element_annotations', 
                      'histone_stemloop_elements.bed.sl3')
    return fp

def hitsclip_discardUnclipped_fq_dir():
    p = os.path.join(get_data_dir(), 
                     'HITS-CLIP_SLBP_histone_mRNA_01', 
                     'adapter_clipped_reads_fastq')
    return p

def hitsclip_bam_dir():
    p = os.path.join(get_data_dir(), 
                     'HITS-CLIP_SLBP_histone_mRNA_01', 
                     'bwa_samse_hg19')
    return p

def hitsclip_fq_fp_l():
    # getting relevant sample data filepaths
    fp_l = sorted(os.listdir(hitsclip_discardUnclipped_fq_dir()))
    fp_l = [os.path.join(hitsclip_discardUnclipped_fq_dir(), fp) for fp in fp_l]
    fq_fp_l = [fp for fp in fp_l if fp.split('.')[-2:] == ['fq','gz']]
    print([os.path.basename(s) for s in fq_fp_l]) #debugging
    ##NOTE: expected sort order is:
    #['SLBP_CLIP_high_MW_high_MNase.HISTONLY.discardUnclipped.fq.gz',
    # 'SLBP_CLIP_high_MW_low_MNase.HISTONLY.discardUnclipped.fq.gz', 
    # 'SLBP_CLIP_low_MW_high_MNase.HISTONLY.discardUnclipped.fq.gz', 
    # 'SLBP_CLIP_low_MW_low_MNase.HISTONLY.discardUnclipped.fq.gz', 
    # 'mock_CLIP_high_MW_high_MNase.HISTONLY.discardUnclipped.fq.gz', 
    # 'mock_CLIP_low_MW_high_MNase.HISTONLY.discardUnclipped.fq.gz']
    return fq_fp_l

def hitsclip_cleav_db_l():
    #NOTE: run test_build_readid_db must generate readid databases
    fp_l = sorted(os.listdir(hitsclip_discardUnclipped_fq_dir()))
    fp_l = [os.path.join(hitsclip_discardUnclipped_fq_dir(), fp) for fp in fp_l]
    cleav_db_l = [fp for fp in fp_l if fp.split('.')[-1] == 'readids']
    print([os.path.basename(s) for s in cleav_db_l]) #debugging
    ##NOTE: expected sort order is:
    #['SLBP_CLIP_high_MW_high_MNase.HISTONLY.discardUnclipped.fq.readids',
    # 'SLBP_CLIP_high_MW_low_MNase.HISTONLY.discardUnclipped.fq.readids', 
    # 'SLBP_CLIP_low_MW_high_MNase.HISTONLY.discardUnclipped.fq.readids', 
    # 'SLBP_CLIP_low_MW_low_MNase.HISTONLY.discardUnclipped.fq.readids', 
    # 'mock_CLIP_high_MW_high_MNase.HISTONLY.discardUnclipped.fq.readids', 
    # 'mock_CLIP_low_MW_high_MNase.HISTONLY.discardUnclipped.fq.readids']
    return cleav_db_l

def hitsclip_bam_fp_l_discardUnclipped():
    fp_l = sorted(os.listdir(hitsclip_bam_dir()))
    fp_l = [os.path.join(hitsclip_bam_dir(), fp) for fp in fp_l]
    bam_fp_l = [fp for fp in fp_l if fp.split('.')[-1] == 'bam']
    bam_fp_l_discardUnclipped = [fp for fp in bam_fp_l if fp.split('.')[2] == 'discardUnclipped']
    print([os.path.basename(s) for s in bam_fp_l_discardUnclipped]) #debugging
    ##NOTE: expected sort order is:
    #['SLBP_CLIP_high_MW_high_MNase.HISTONLY.discardUnclipped.bam',
    # 'SLBP_CLIP_high_MW_low_MNase.HISTONLY.discardUnclipped.bam', 
    # 'SLBP_CLIP_low_MW_high_MNase.HISTONLY.discardUnclipped.bam', 
    # 'SLBP_CLIP_low_MW_low_MNase.HISTONLY.discardUnclipped.bam', 
    # 'mock_CLIP_high_MW_high_MNase.HISTONLY.discardUnclipped.bam', 
    # 'mock_CLIP_low_MW_high_MNase.HISTONLY.discardUnclipped.bam']
    return bam_fp_l_discardUnclipped

def hitsclip_bam_fp_l():
    fp_l = sorted(os.listdir(hitsclip_bam_dir()))
    fp_l = [os.path.join(hitsclip_bam_dir(), fp) for fp in fp_l]
    bam_fp_l = [fp for fp in fp_l if fp.split('.')[-1] == 'bam']
    bam_fp_l = [fp for fp in bam_fp_l if fp.split('.')[2] != 'discardUnclipped']
    print([os.path.basename(s) for s in bam_fp_l]) #debugging
    ##NOTE: expected sort order is:
    #['SLBP_CLIP_high_MW_high_MNase.HISTONLY.bam',
    # 'SLBP_CLIP_high_MW_low_MNase.HISTONLY.bam',
    # 'SLBP_CLIP_low_MW_high_MNase.HISTONLY.bam',
    # 'SLBP_CLIP_low_MW_low_MNase.HISTONLY.bam',
    # 'mock_CLIP_high_MW_high_MNase.HISTONLY.bam',
    # 'mock_CLIP_low_MW_high_MNase.HISTONLY.bam']
    return bam_fp_l

def hitsclip_n_mapped_reads_l():
    # the normalization factor is parameterized because user sometimes 
    # wants to evaluate a subsetted bam file but also desires normalization by
    # the total number of mapped reads from the upstream bam file
    # e.g. the HISTONLY files are subsetted and do not contain all reads
    n_mapped_reads_l = [17931502, #SLBP_CLIP_high_MW_high_MNase.bam
                        14698879, #SLBP_CLIP_high_MW_low_MNase.bam
                        20020463, #SLBP_CLIP_low_MW_high_MNase.bam
                        15256431, #SLBP_CLIP_low_MW_low_MNase.bam
                        18992515, #mock_CLIP_high_MW_high_MNase.bam
                        21083226, #mock_CLIP_low_MW_high_MNase.bam
                       ]
    return n_mapped_reads_l

def remove_test_files():
    # Here are lists of filenames generated by the tests. This is essentially
    # the same files as listed in the git ignore file.
    bg_fn_l = os.listdir(hitsclip_bam_dir())
    bg_fn_l = [s for s in bg_fn_l if s.split('.')[-1] == 'bg']
    bg_fn_l.sort()
    print('Removing bedgraph files that were generated during tests...')
    for fn in bg_fn_l:
        fp = os.path.join(hitsclip_bam_dir(), fn)
        if os.path.isfile(fp):
            print(fp)
            os.remove(fp)
    
    readids_db_fn_l = os.listdir(hitsclip_discardUnclipped_fq_dir())
    readids_db_fn_l = [s for s in readids_db_fn_l if s.split('.')[-1] == 'readids']
    readids_db_fn_l.sort()
    print('Removing readid db files that were generated during tests...')
    for fn in readids_db_fn_l:
        fp = os.path.join(hitsclip_discardUnclipped_fq_dir(), fn)
        if os.path.isfile(fp):
            print(fp)
            os.remove(fp)
    
    print('Removing sl3 cis-element file that was generated during tests...')
    fp = histone_sl_bed_sl3_fp()
    if os.path.isfile(fp):
        print(fp)
        os.remove(fp)
    
    return

