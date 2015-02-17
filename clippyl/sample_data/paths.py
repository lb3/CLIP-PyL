"""This module contains numerous function that produce 
relative paths to the sample data"""
import os

def get_data_dir():
    data_dir = os.path.dirname(os.path.abspath(__file__))
    return data_dir

def histone_gene_bed_fp():
    return os.path.join(get_data_dir(), 'genomic_interval_queries', 'histone_genes_HGNC.bed')

def histone_sl_bed_fp():
    return os.path.join(get_data_dir(), 'cis_element_annotations', 'histone_SL_elements.bed')

def histone_sl_bed_sl3_fp():
    return os.path.join(get_data_dir(), 'cis_element_annotations', 'histone_SL_elements.bed.sl3')

def hitsclip_discardUnclipped_fq_dir():
    return os.path.join(get_data_dir(), 'HITS-CLIP_SLBP_histone_mRNA_01', 'adapter_clipped_reads_fastq')

def hitsclip_bam_dir():
    return os.path.join(get_data_dir(), 'HITS-CLIP_SLBP_histone_mRNA_01', 'bwa_samse_hg19')

#TODO: make a clean-up test files command for files created by tests
