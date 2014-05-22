import os

def histone_gene_bed_fp():
    return os.path.join(os.path.getcwd(), 'histone_genes_HGNC.bed')

def histone_sl_bed_fp():
    return os.path.join(os.path.getcwd(), 'histone_SL_elements.bed')


