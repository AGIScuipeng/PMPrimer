# PMPrimer - Automatically design multiplex PCR primer pairs for diverse templates
中文说明, [请点击](README_CN.md)
## 0x0 Module Information
> This software need third package : primer3-py, biopython, seaborn etc.

> If want use BLAST, need local support ncbi-blast-2.13.0+.

> This software running on Linux, not mean it cannot work on Windows.
+ Input file
+ Basic process for data
+ Primer design and extract of multiple alignment
+ Amplicon design and evaluate
+ Debug mode

## 0x1 Installation
1. Install by code
   
   Download source code from this repo, and install primer3-py, biopython, seaborn in your devices, then you can run this software by use `python3 pmprimer.py`

2. Install by pip
   
   Install this software by use command `pip(or pip3 for python3) install PMPrimer`.
   
   If u use pip to install this software, pip will automatically install dependency packages and add a command pmprimer in your system.

## 0x2 Usage
1. Input file
   > --file -f

    Use like : `pmprimer -f seqs.fasta or --file seqs.fasta`

2. Basic process for data
   > --progress -p

    All parameter as follows :
    + notlen Do not processing sequences according to distribution of length
    + notsameseq Do not remove redanduncy sequences according to three levels
    + matrix Draw heatmap according to mismatch analysis.

    >If length distribution is good, should use 'notlen' prevent removal of minority sequences.

    Add parameters behind -p, use like : `-p notlen notsameseq matrix`, order doesn't matter.

3. Primer design and extract of multiple alignment
   >--alldesign -a

    All parameter as follows :
    + threshold Threshold for identify conservative region, default is 0.95, use like : threshold:0.95
    + minlen Minimum continues length when identify conservative regions, default is 15bp, use like : minlen:15
    + gaps Maximum rate of gaps, default is 0.1, use like : gaps:1.0
    + tm Default TM, default is 50, use like : tm:50.0
    + muscle Use MUSCLE to align sequences and save as another file
    + merge Use merge module to merge conservative regions
    + primer2 Primer design and extract
    + pdetail2 Information of primer2

    add parameters behind -a, like : `-a threshold:0.96 minlen:16 merge primer2`, order doesn't matter.

4. Amplicon design and evaluate
   > --evaluate -e

    All parameter as follows :
    + hpcnt Maximum count of haploType, default is 10, use like : hpcnt:10
    + minlen Minimum length of amplicon, default is 150bp, use like : minlen:150
    + maxlen Maximum length of amplicon, default is 1500bp, use like : maxlen:1500
    + blast Use fasta to create blast db and search, use ',' split different files, use like : blast:../taxo1.fasta,homo.fasta
    + rmlow Remove TM lower than configure in Design module
    + save Save final results

    add parameters behind -e, like : `-e hpcnt:50 maxlen:1200 save`, order doesn't matter.

5. Debug mode
   > --debuglevel -d

    parameter is number :
    + 1 Basic debug info
    + 2 Mocule debug info
    + 4 Complex debug info

    If need debug info, use `-d 1` in normal conditions, this number use binary calculate.

## 0x3 Tips
1. If input file is alignments, use like :

    `pmprimer -f seqs.fasta -a merge primer2 -e hpcnt:70 save`

2. If input file need data process, use like :

    `pmprimer -f seqs_need_filt.fasta -p default`
    > this command will process data according to major length and different third level when same sequences.

3. If input file need align, use like : 

    `pmprimer -f seqs_need_align.fasta -a muscle merge primer2 -e hpcnt:70 save`
    > this command will save alignments as 'seqs_need_align.mc.fasta'

## 0x4 Datasets In Paper

Dataset in paper can obtained by https://github.com/AGIScuipeng/PMPrimer_datasets

1. 16S ribosomal RNA (rRNA) genes of Archaea
   
   Archaea_16SrRNA.rep.mc.fasta is more than 200MB, so only upload original file to github, need use use MUSCLE5 to align sequences and save output file as `Archaea_16SrRNA.rep.mc.fasta` or use `pmprimer -f Archaea_16SrRNA.rep.fasta -p notlen notsameseq muscle` to align sequences.
   
   Command in paper is : `pmprimer -f Archaea_16SrRNA.rep.mc.fasta -a threshold:0.85 gaps:1.0 merge primer2 haplo tm:45.0 -e hpcnt:600 save`

> 
2. hsp65 (groEL2) genes of Mycobacteriaceae
   
   Command in paper is : `pmprimer -f Mycobacteriaceae_groEL2.filt.mc.fasta -a primer2 -e hpcnt:70 save`
> 
3. tuf genes of Staphylococci
   
   Command in paper is : `pmprimer -f Staphylococcus_tuf.filt.mc.fasta -a threshold:0.995 minlen:5 merge primer2 -e save`