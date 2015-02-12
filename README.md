#CLIP-PyL

This is the ALPHA release of the CLIP-PyL package, which contains scripts that can parse "crosslink signatures" from aligned CLIP-seq data.

Currently, HITS-CLIP is the only CLIP-seq method that is supported. However, PAR-CLIP and iCLIP data will be supported if considerable interest arises.

This tool requires that the user pre-processes the "raw" HITS-CLIP reads by clipping the adapter sequences. I typically use the FASTQ/A Clipper tool from the [FASTX toolkit](hannonlab.cshl.edu/fastx_toolkit/) for this task. However, other suitable options exist. The user must also align the preprocessed reads to a reference genome assembly. Currently, single-end read alignments generated by the [bwa](bio-bwa.sourceforge.net/bwa.shtml) fast read alignment software are supported. Finally, the user must use the [SAMTools software](http://www.htslib.org/) to generate indexed bam files. CLIP-PyL takes the indexed bam files and returns a set of [bedgraph](http://genome.ucsc.edu/FAQ/FAQformat.html#format1.8) files containing basewise crosslink metrics. The crosslink metrics that are generated include single base-pair deletions, fragment termini and raw coverage. Note that the bedgraph files can be viewed with most genome browser tools (e.g. [Broad's Integrated Genome Viewer](http://www.broadinstitute.org/igv/) or [UCSC's Genome Browser](http://genome.ucsc.edu/)).

The CLIP-PyL package also contains utilities to produce pdf files containing coverage map graphics for sets of genome intervals, which are specified by the user in a [bed file format](http://genome.ucsc.edu/FAQ/FAQformat.html#format1). Notably, cis-element regions can be included in the coverage graphics if the user supplies additional bed files specifying the locations of these elements. CLIP-PyL utilizies the third-party graphical back-end known as [matplotlib](http://matplotlib.org/) to generate these coverage graphics.

#Dependencies

* python >= 3.2 https://www.python.org/

* pysam https://pypi.python.org/pypi/pysam

* matplotlib http://matplotlib.org/

If you have multiple python versions installed on your machine then please be certain to install the pysam and matplotlib packages into your python3 environment.

#Installation

Download the [current release version of the package]().

The clippyl package has a command line interface (CLI). The clippyl CLI can be invoked by command liscripts are contained within the clippyl package. The you can test your installation by issuing the following command:

python3 clippyl

#Using CLIP-PyL

    #!/usr/bin/env python3
    
    from towelstuff import location
    from towelstuff import utils
    
    if utils.has_towel():
        print "Your towel is located:", location.where_is_my_towel()



