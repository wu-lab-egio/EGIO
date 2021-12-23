## EGIO

EGIO (Exon Group Ideogram based detection of Orthologous exons and Orthologous isoforms) is aimed to detect orthologous exons and isoforms within pre-defined orthologous gene pairs. 

## requirement
(1) EGIO required python packages pandas and numpy, to install pandas and numpy in terminal

  pip install pandas
  pip install pandas

(2) To increase the accuracy, the algorithm uses a BLASTN guided model, so a reciprocal BLASTN of exons is required. To run recipral BLASTN, BLAST+ is required before running EGIO, which could be found in https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/. An example is listed here:

    hsa	ptr
  
    ENSG00000131018	ENSPTRG00000018717
  
    ENSG00000183091	ENSPTRG00000012536
  
    ENSG00000183091	ENSPTRG00000012536
  
  ...

  
(3) a tab seperated file containind pre-defined orthologous gene pairs, which could be prepared according to Inparanoid: https://inparanoid.sbc.su.se/cgi-bin/index.cgi. The gene id is Uniprot ID, which can be transformed to Ensembl Gene ID using gene annotations in Ensembl BioMart.

(4) reference files, which can be downloaded from Ensembl (http://asia.ensembl.org/info/data/ftp/index.html):
  cDNA fasta files,
  CDS fasta files,
  reference transcriptome gtf files
  
  
## Usage

