import time
import argparse
import sys
import os

from clipPyL.flatfile_parsing import FastqReader, validate_fastq_format
from clipPyL.chrom_db_io import HitsClipSQLite

def build_db(in_sam_fp, in_cleav_fp, out_dir, clipped_read_filetype='fastq'):
    
    in_sam_fn = os.path.basename(in_sam_fp)
    in_sam_fn_prefix = os.path.splitext(in_sam_fn)[0]
    
    #TODO: give option to read a list of read Ids instead of fastq file (write new
    #      input_readID_list method in HitsClipSQLite class)
    #TODO: parameterize uniq_only in argparse below...note this may be overridden 
    #      by permute routine, as all reads must be read into the sam file to sample
    #      from the full population of alignments. give a permute paramenter that
    #      can be None, all, unique
    #TODO: save temporary files to a subdirectory and delete after you're done
    #      give option to clean-up files after the script has run, with default as
    #      yes to cleanup
    #TODO: detect file if it already exists and warn of overwrite (write a 
    #      function to call at each HitsClipSQLite file creation)
    
    if in_cleav_fp == in_sam_fp:
        print('#######################################')
        print('''You indicated that the alignment file is comprised''')
        print('''solely of adapter-clipped reads.''')
        print('#######################################')
        in_cleav_fp = None 
        #setting to None triggers proper build_basewise_signature routine
        #TODO: populate n_of_cleaved_reads
    
    elif clipped_read_filetype == 'fastq':
        
        print('#######################################')
        in_cleav_bn = os.path.basename(in_cleav_fp)
        cleav_db_fp = os.path.join(out_dir, in_cleav_bn + '.sqlite')
        
        print('Validating fastq...')
        if validate_fastq_format(in_cleav_fp):
            pass
        else:
            print('something bad happened')
            raise IOError
        
        cleav_db = HitsClipSQLite(cleav_db_fp)
        
        start_time = time.time()
        n_of_cleaved_reads = cleav_db.input_fastq(in_cleav_fp)
        cleav_db.conn.close()
        elapsed_time = time.time() - start_time
        
        print(n_of_cleaved_reads)
        print('Adapter-clipped read IDs were written and indexed to:')
        print(cleav_db_fp)
        print('The amount of time that elapsed during the process was:')
        print('{0:.2f}'.format(round(elapsed_time,2)) + ' seconds')
        print('#######################################')
    
    elif clipped_read_filetype == 'list':
        #TODO: implement me and write a function in fastq_reader to
        #      generate this list from a fastq file
        #TODO: populate n_of_cleaved_reads
        raise IOError
    
    else:
        print('something bad happened')
        raise IOError
    
    print('#######################################')
    sam_db_fp = os.path.join(out_dir, in_sam_fn + '.sqlite')
    
    start_time = time.time()
    sam_db = HitsClipSQLite(sam_db_fp)
    sam_stat_dict = sam_db.input_sam(in_sam_fp)
    # note: uniqs are selected downstream during basewise signature 
    # calculation (see below) so that all reads are present in the db 
    # for permutation of the 1D offsets
    sam_db.conn.close()
    
    elapsed_time = time.time() - start_time
    print('bwa-samse alignments were written and indexed to:')
    print(sam_db_fp)
    print('The amount of time that elapsed during the process was:')
    print('{0:.2f}'.format(round(elapsed_time,2)) + ' seconds')
    print('#######################################')
    
    print('#######################################')
    start_time = time.time()
    
    pyl_db_fp = os.path.join(out_dir, in_sam_fn_prefix + '.clipPyL.sqlite')
    #pyl_db_fp is a sqlite datbase that holds the basewise clip signature
    
    #TODO: parametrize uniq_only in command line interface
    hitsclip_pyl_SQLite = HitsClipSQLite(pyl_db_fp)
    bw_stat_dict = hitsclip_pyl_SQLite.build_basewise_signature( 
                                          sam_db_fp = sam_db_fp,
                                          cleav_db_fp = cleav_db_fp, 
                                          uniq_only = True )
    
    hitsclip_pyl_SQLite.conn.close()
    
    elapsed_time = time.time() - start_time
    print('positional hits-clip signature was written and indexed to:')
    print(pyl_db_fp)
    print('The amount of time that elapsed during the process was:')
    print('{0:.2f}'.format(round(elapsed_time,2)) + ' seconds')
    print('#######################################')
    
    #TODO: SAVE THESE STATS TO A TABLE IN THE DB FOR THE PURPOSES OF NORMALIZATION
    # WITH DOWNSTREAM METHODS
    stat_fp = os.path.join(out_dir, in_sam_fn_prefix + '.clipPyL.log')
    with open(stat_fp, 'w') as stat_fh:
        stat_fh.write('hits-clip PyL was built with the following command line parameters:\n')
        stat_fh.write('--alignments: ' + in_sam_fp + '\n')
        stat_fh.write('--clipped_reads: ' + in_cleav_fp + '\n')
        stat_fh.write('--out_dir: ' + out_dir + '\n')
        stat_fh.write('--clipped_read_filetype: ' + str(clipped_read_filetype) + '\n')
        stat_fh.write('\n')
        stat_fh.write('File stats:\n')
        stat_fh.write('Total number of adapter-clipped reads: ' + str(n_of_cleaved_reads) + '\n')
        #total number of reads is gleaned from sam file
        stat_fh.write('Total number of reads: ' + str(sam_stat_dict['n_of_reads']) + '\n')
        #gleaned from sam flagbits
        stat_fh.write('Number of mapped reads: ' + str(sam_stat_dict['n_of_mapped_reads']) + '\n')
        #gleaned from sam XT opt field
        stat_fh.write('Number of uniquely mapped reads: ' + str(sam_stat_dict['n_of_unique_aligns']) + '\n')
        # total of nt sites with coverage
        stat_fh.write('Number of nucleotides with coverage: ' + str(bw_stat_dict['n_of_nts_covered']) + '\n')
        # total number of nts with oneD
        stat_fh.write('Number of nucleotides with 1D: ' + str(bw_stat_dict['n_of_oneD_sites']) + '\n')
        #total number of nts with cleavage
        stat_fh.write('Number of nucelotides with cleavage: ' + str(bw_stat_dict['n_of_cleavage_sites']) + '\n')
        #note: impacted by uniq_only
        stat_fh.write('Number of mapped adapter-clipped reads: ' + str(bw_stat_dict['n_of_mapped_cleav_reads']) + '\n')
        #note: impacted by uniq_only
        stat_fh.write('Number of mapped 1D reads: ' + str(bw_stat_dict['n_of_mapped_oneD_reads']) + '\n')
    
    return


class Usage(Exception):
    def __init__(self, exitStat):
        self.exitStat = exitStat

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            d = '''Calculate the base-wise HITS-CLIP signature and buld the '.PyL.sqlite' database'''
            parser = argparse.ArgumentParser(description = d)
            
            parser.add_argument('-a', '--alignments', 
                                metavar='file path to alignments', type=str,
                                required=True,
                                help='''Provide the filepath to the read alignments from your HITS-CLIP
                                        experiment. This must point to a sam file generated by the 
                                        bwa samse algorithm.''')
            
            parser.add_argument('-c', '--clipped_reads', 
                                metavar='file path to adapter-clipped reads', type=str,
                                required=True,
                                help='''Provide the filepath to the read IDs of the adapter-clipped 
                                        read subset. The expected file format is fastq.\n
                                        However, if you enter the same filepath that you entered for
                                        the -a (--alignemnts) argument then it is assumed that the 
                                        alignment file is comprised solely of adapter-clipped reads.\n
                                        Also, see --in_cleav_ft argument for another option.''')
            
            parser.add_argument('-o', '--output_dir', 
                                metavar='directory path for output', type=str,
                                required=False, default=os.getcwd(),
                                help='''Provide an output directory. Default output directory is
                                        the current working directory''')
            
            parser.add_argument('--clipped_read_filetype', 
                                metavar='format type of clipped read file', type=str, 
                                required=False, default='fastq', choices=['fastq','list'],
                                help='''Set this parameter to list if your adpater-clipped read file is 
                                        a list of read IDs. By default, a fastq file is expected''')
            
            args = parser.parse_args(argv[1:])
            
        except SystemExit as exitStat:
            raise Usage(exitStat)
        
        print(args.alignments)
        print(args.clipped_reads)
        print(args.output_dir)
        print(args.clipped_read_filetype)
        
        #TODO: WHY DO THESE args return a list object?! try to constrain to return a string
        build_db(in_sam_fp=args.alignments, 
                 in_cleav_fp=args.clipped_reads, 
                 out_dir=args.output_dir, 
                 clipped_read_filetype=args.clipped_read_filetype)
    
    except Usage as err:
        return err.exitStat

