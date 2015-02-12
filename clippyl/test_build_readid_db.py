import os

from clippyl.build_readid_db import build_ReadidSQLite_dbs
from clippyl.sample_data.paths import hitsclip_discardUnclipped_fq_dir

if __name__ == '__main__':
    # getting relevant sample data filepaths
    fp_l = os.listdir(hitsclip_discardUnclipped_fq_dir())
    fp_l = [os.path.join(hitsclip_discardUnclipped_fq_dir(), fp) for fp in fp_l]
    fq_fp_l = [fp for fp in fp_l if fp.split('.')[-2:] == ['fq','gz']]
    print (fq_fp_l)
    
    build_ReadidSQLite_dbs(fq_fp_l)
