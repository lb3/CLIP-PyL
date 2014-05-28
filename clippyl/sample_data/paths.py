import os

def get_data_dir():
    data_dir = os.path.dirname(os.path.abspath(__file__))
    return data_dir

def histone_gene_bed_fp():
    return os.path.join(get_data_dir(), 'histone_genes_HGNC.bed')

def histone_sl_bed_fp():
    return os.path.join(get_data_dir(), 'histone_SL_elements.bed')


