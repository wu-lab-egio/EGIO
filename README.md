# Introduction

EGIO (Exon Group Ideogram based detection of Orthologous exons and Orthologous isoforms) is aimed to detect orthologous exons and isoforms within pre-defined orthologous gene pairs. EGIO uses a dynamic programming strategy to do isoform alignment, with a reciprocal BLASTN result to guide the whole process.

# Requirement
(1) EGIO required python packages pandas and numpy, to install pandas and numpy in terminal

    pip install pandas
  
    pip install numpy


(2) To increase the accuracy, the algorithm uses a BLASTN guided model, so a reciprocal BLASTN of exons is required. To run recipral BLASTN, BLAST+ is required before running EGIO, which could be found in https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/. An example is listed here:

    hsa	ptr
  
    ENSG00000131018	ENSPTRG00000018717
  
    ENSG00000183091	ENSPTRG00000012536
  
    ENSG00000183091	ENSPTRG00000012536
  
    ... 

To be noted, the header is required in the file.
 

(3) a tab seperated file containind pre-defined orthologous gene pairs, which could be prepared according to Inparanoid: https://inparanoid.sbc.su.se/cgi-bin/index.cgi. The gene id is Uniprot ID, which can be transformed to Ensembl Gene ID using gene annotations in Ensembl BioMart.


(4) reference files, which can be downloaded from Ensembl (http://asia.ensembl.org/info/data/ftp/index.html):

  cDNA fasta files,
  
  CDS fasta files,
  
  reference transcriptome gtf files
  
  
# Usage
To test EGIO, please download example files in https://github.com/EGIO1/EGIO_example_source, and unzip gtf files. Put these files in folder "example", and put the "example" folder into "EGIO-main". When running EGIO, it may notice that it it not permitted to run _Run_egio.sh. 

To run EGIO, type following commands in Ternimal:


    cd /path/to/EGIO-main
    chmod 777 _Run_egio.sh
    ./_RUN_egio.sh -s hsa -S ptr -r example/hsa.gtf -R example/ptr.gtf -e example/hsa_mRNA_example.fa -E example/ptr_mRNA_example.fa -o example/hsa_CDS_example.fa -O example/ptr_CDS_example.fa -h example/homogene.txt -p 6 -i 0.8 -c 0.8 -m 2 -n -2 -g -1

For macOS, it might also notice that it it not permitted to run ___pairwisealign.so, 
    
    chmod 777 ___pairwisealign.so



Explaination:

-s: name of species1, to be noted, the name should be consistent with the header of files containing pre-defined orthologous gene pairs, eg: hsa (required)

-S: name of species2, to be noted, the name should be consistent with the header of files containing pre-defined orthologous gene pairs, eg: ptr (required)

-r: gtf file of species1 (required)

-R: gtf file of species2 (required)

-e: cDNA fasta file of species1 (required)

-E: cDNA fasta file of species2 (required)

-o: CDS fasta file of species1 (required)

-O: CDS fasta file of species2 (required)

-h: orthologous/homologous gene pair containing file (required)


-p: multiple processing to run EGIO (optional, default is 1)

-i: identity threshold to detect orthologous exons during dynamic programming (optional, default is 0.8)

-c: coverage threshold to detect orthologous exons during dynamic programming (optional, default is 0.8)

-m: match score during dynamic programming (optional, default is 2)

-n: mismatch penalty during dynamic programming (optional, default is -2)

-g: gap penalty during dynamic programming (optional, default is -1)


# Description of EGIO results

Two files will be generated after running EGIO: ExonGroup.txt and OrthoIso.txt.

## Explaination of ExonGroup.txt headers:

Group: exon group number

species1EnsemblG: gene id of species1

species1Pos: unique exon region of a exon group in species1

species2EnsemblG: gene id of species2

species2Pos: unique exon region of a exon group in species2

Iden: identity of unique exon region in two species of a exon group

Type: type of exon group.


## Explaination of OrthoIso.txt headers:

species1: gene id of species1

species2: gene id of species1

species1iso: isoform of species1

species2iso: isoform of species2

exoniden: identity of corresponding exons
