# PMPrimer: a Python-based tool for automated design and evaluation of multiplex PCR primer pairs for diverse templates

Unlike single sequence, primer design on multiple sequences needs data quality filtering, multiple sequence alignment, conservative region identification, primer design for multiple templates, and evaluation of coverage and taxon specificity on multiple templates.

PMPrimer is a Python-based tool for automated design and evaluation of multiplex PCR primer pairs for diverse templates. The greatest strength of PMPrimer is the ability to identify conservative region using Shannon's entropy method, tolerate gaps using haplotype-based method, and evaluate multiplex PCR primer pairs according to template coverage and taxon specificity. PMPrimer were tested in datasets with diversified conservation and data size, including *tuf* genes of Staphylococci, *hsp65* genes of Mycobacteriaceae, and 16S ribosomal RNA genes of Archaea. By comparison with previously designed primers and existing tools, PMPrimer has outstanding performance, such as automation, big data support, higher accuracy, and reasonable running time. 

PMPrimer is a command-line tool. To handle diversified characteristics of target sequences, PMPrimer has multiple built-in parameters to make appropriate adjustments for data processing. The only necessary input file for PMPrimer is multiple sequences of target gene in FASTA format. Results are provided as JSON and CSV format files, which can be easily loaded into other programming languages for further analysis.

## 0x01 Installation
The software requires third-party packages as follow: primer3-py, biopython, matplotlib, pandas, seaborn and blast (2.13.0).

All python-based third-party packages will be automatically installed by pip. User needs install blast according to the documentation at [BLAST OFFICIAL WEBSITE](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs).

PMPrimer has been tested with Python3.8+ on Linux and Windows.

1. Install by pip

    We recommend that you install PMPrimer through pip using the command: `pip install PMPrimer` or `pip3 install PMPrimer`. pip will automatically install all dependency packages and add a command "pmprimer" in your system.

2. Install by code

    Download PMPrimer source code from this repository, decompress the source code, install all dependency packages by yourself in your system, then you can run this software by use `python3 PMPrimer/pmprimer.py`.

## 0x02 Usage

Use help command `pmprimer --help` to show usage message. The detailed information is provided as follow.

1. Input file
   > --file -f

    The only necessary input file for PMPrimer is multiple sequences of target gene in FASTA format. It can be unaligned or aligned. If you assign an unaligned file, you need assignment `-a muscle` parameter in primer design model to do multiple sequence alignment.

    Use like : `pmprimer -f seqs.fasta` or `pmprimer --file seqs.fasta`

2. Basic process
   > --progress -p

    In basic process, PMPrimer assesses data quality, calculates length distribution of input data, filters low quality templates (such as too small and too long templates) according to length distribution, removes redundant templates in terminal taxa in default. If input file has been manually checked, you can use `-p notlen` to disable length filter based on length distribution. Similarly, you can use `-p notsameseq` to disable redundant filter based on sequence comparison in terminal taxa. Moreover, if you want to draw heatmap according to distances of multiple templates, you can use `-p matrix` to trigger drawing heatmap according to distances of multiple templates.

    All parameter as follows :
    + notlen Do not processing sequences according to distribution of length
    + notsameseq Do not remove redundancy sequences in terminal taxa
    + matrix Draw heatmap according to distances of multiple templates

    Add parameters behind -p, use like : `-p notlen notsameseq matrix`, order doesn't matter.

3. Primer design
   >--alldesign -a

    This module includes three steps, multiple sequence alignment by `muscle` parameter, conservative region identification by `threshold minlen gaps merge` parameters, and primer design by `tm primer2` parameters. If you assign an unaligned file, you must asign `muscle` parameter in multiple sequence alignment step. You can adjust parameters for conservative region identification to obtain optimal conservative regions. You need assign `primer2` to trigger primer design when you get optimal conservative regions.

    All parameter as follows :
    + muscle Use MUSCLE to align multiple sequences, the alignment sequences will be saved as *.mc.fasta
    + threshold Threshold for identify conservative region, default is 0.95 frequency for major allele, use like threshold:0.95
    + minlen Minimum length for identifying conservative region, default is 15, use like minlen:15
    + gaps Maximum rate of gaps, default is 0.1, use like gaps:0.1
    + merge Use merge module to merge conservative regions
    + rank1 Display the diversity score ranking of non conservative regions
    + rank2 Display the diversity score ranking of conservative regions
    + haplo Display the count of haploType sequences of conservative regions
    + tm Minimum melting temperature, default is 50, use like tm:50.0
    + primer2 Primer design

    Add parameters behind -a, use like : `-a threshold:0.96 minlen:16 merge primer2`, order doesn't matter.

4. Amplicon selection and evaluation
   > --evaluate -e

    To evaluate amplicon specificity, we use blast to align amplicon primer pairs to target fasta files. We can assign multiple target fasta files by using "," to split file pathes.

    All parameter as follows :
    + hpcnt Maximum count of haploType sequences, default is 10, use like hpcnt:10
    + minlen Minimum length of amplicon, default is 150bp, use like : minlen:150
    + maxlen Maximum length of amplicon, default is 1500bp, use like : maxlen:1500
    + blast Use blast to evaluate amplicon specificity by aligning amplicon primer pairs to target fasta files, use ',' to split multiple target files, use like blast:seqs.fasta,hosts.fasta
    + rmlow Remove primers with melting temperature lower than parameter "tm" in primer design module
    + save Save final results

    Add parameters behind -e, use like : `-e hpcnt:50 maxlen:1200 save`, order doesn't matter.

5. Debug mode
   > --debuglevel -d

    All parameter is number :
    + 1 Basic debug info
    + 2 Module debug info
    + 4 Complex debug info

    If need debug info, use like : `-d 1` in normal conditions, this number use binary calculate.

## 0x03 Tips
1. If input file need data process, use like :

    `pmprimer -f seqs.fasta -p default`
    > Filt sequences according to distribution of length,remove redundancy sequences in terminal taxa.

2. If input file is alignments, use like :

    `pmprimer -f seqs.fasta -a primer2`

3. If input file need align, use like : 

    `pmprimer -f seqs.fasta -a muscle merge primer2`
    > This command will save alignments as 'seqs.mc.fasta'

4. You can use `rank1`, `rank2`, `haplo` parameters to check if the conservative region meets your requirements
   
    `pmprimer -f seqs.fasta -a rank1 rank2 haplo`
    
    If you think the current conservative regions is optimal, use `primer2` to trigger primer design

## 0x04 Datasets In Paper

Dataset in paper can obtained from [PMPrimer Datasets](https://github.com/AGIScuipeng/PMPrimer_datasets)

1. 16S ribosomal RNA (rRNA) genes of Archaea
   
   Command in paper is : `pmprimer -f Archaea_16SrRNA.rep.mc.fasta -a threshold:0.85 gaps:1.0 merge primer2 haplo tm:45.0 -e hpcnt:600 save`

> 
1. hsp65 (groEL2) genes of Mycobacteriaceae
   
   Command in paper is : `pmprimer -f Mycobacteriaceae_groEL2.filt.mc.fasta -a primer2 -e hpcnt:70 save`
> 
1. tuf genes of Staphylococci
   
   Command in paper is : `pmprimer -f Staphylococcus_tuf.filt.mc.fasta -a threshold:0.995 minlen:5 merge primer2 -e save`