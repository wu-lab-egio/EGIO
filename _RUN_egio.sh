#!/bin/bash

# EGIO pipline is based on python scripts and requirs following programs and python packages
# 1. gcc
# 2. BLAST+ program [https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/]
# 3. python packages: pandas, numpy

while getopts ":s:S:r:R:e:E:o:O:h:i:c:m:n:g:p:" opt; do
  case $opt in
    s) species1="$OPTARG"
    ;;
    S) species2="$OPTARG"
    ;;
    r) GTFreference1="$OPTARG"
    ;;
    R) GTFreference2="$OPTARG"
    ;;
    e) cdna1="$OPTARG"
    ;;
    E) cdna2="$OPTARG"
    ;;
    o) cds1="$OPTARG"
    ;;
    O) cds2="$OPTARG"
    ;;
    h) ortholog="$OPTARG"
    ;;
    i) identity="$OPTARG"
    ;;
    c) coverage="$OPTARG"
    ;;
    m) match="$OPTARG"
    ;;
    n) mismatch="$OPTARG"
    ;;
    g) gap="$OPTARG"
    ;;
    p) pnum="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

## setting parameters if missing in command
if [ ! -z "${species1}" ];then
  echo "species1 is $species1"
else
  echo "species1 is not set, which will be set as egiospe1"; 
  species1=egiospe1
fi

if [ ! -z "${species2}" ];then
  echo "species2 is $species2"
else
  echo "species2 is not set, which will be set as egiospe2"; 
  species2=egiospe2
fi

if [ -z "${GTFreference1}" ];then
  echo "GTFreference1 is required"
  exit 1
fi

if [ -z "${GTFreference2}" ];then
  echo "GTFreference2 is required"
  exit 1
fi

if [ -z "${cdna1}" ];then
  echo "cdna1 is required"
  exit 1
fi

if [ -z "${cdna2}" ];then
  echo "cdna2 is required"
  exit 1
fi

if [ -z "${cds1}" ];then
  echo "cds1 is required"
  exit 1
fi

if [ -z "${cds2}" ];then
  echo "cds2 is required"
  exit 1
fi

if [ -z "${ortholog}" ];then
  echo "ortholog is required"
  exit 1
fi

if [ ! -z "${identity}" ];then
  echo "identity is $identity"
else
  echo "identity is required"
  exit 1
fi

if [ ! -z "${coverage}" ];then
  echo "coverage is $coverage"
else
  echo "coverage is required";
  exit 1
fi

if [ ! -z "${match}" ];then
  echo "match score is $match"
else
  ((match = 2))
  echo "match score is $match"
fi

if [ ! -z "${mismatch}" ];then
  echo "mismatch penalty is $mismatch"
else
  ((mismatch = -2))
  echo "mismatch penalty is $mismatch"
fi

if [ ! -z "${gap}" ];then
  echo "gap penalty is $gap"
else
  ((gap = -1))
  echo "gap penalty is $gap"
fi

if [ ! -z "${pnum}" ];then
  echo "use $pnum core"
else
  ((pnum = 1))
  echo "use $pnum core"
fi

###=======================================
if [ ! -f "___pairwisealign.so" ]; then
echo "compile ___pairwisealign.c"
gcc  -fPIC -shared -m64 -o ___pairwisealign.so ___pairwisealign.c
fi
###=======================================
echo "prepare required files of egio"

if [ ! -d "extrainfo" ]; then
mkdir "extrainfo"
fi

if [ ! -d "extrainfo/$species1" ]; then
mkdir "extrainfo/$species1"
fi

if [ ! -f "extrainfo/$species1.exon" ]; then
python __prepare_egio_extra.py --gtf $GTFreference1 --cdna $cdna1 --cds $cds1 --species $species1
fi

if [ ! -f "extrainfo/$species1/$species1.ndb" ]; then
makeblastdb -in extrainfo/$species1.exonfasta -dbtype nucl -out extrainfo/$species1/$species1
fi

###--------------------------------------
if [ ! -d "extrainfo/$species2" ]; then
mkdir "extrainfo/$species2"
fi

if [ ! -f "extrainfo/$species2.exon" ]; then
python __prepare_egio_extra.py --gtf $GTFreference2 --cdna $cdna2 --cds $cds2 --species $species2
fi

if [ ! -f "extrainfo/$species2/$species2.ndb" ]; then
makeblastdb -in extrainfo/$species2.exonfasta -dbtype nucl -out extrainfo/$species2/$species2
fi

###=======================================
if [ ! -f "extrainfo/$species1-$species2" ]; then
echo "reciprocal blastn of $species1 and $species2: $species1 as query"
blastn -query extrainfo/$species1.exonfasta -out extrainfo/$species1-$species2 -db extrainfo/$species2/$species2 -outfmt 6 -evalue 1e-5 -num_threads $pnum
fi

if [ ! -f "extrainfo/$species2-$species1" ]; then
echo "reciprocal blastn of $species1 and $species2: $species2 as query"
blastn -query extrainfo/$species2.exonfasta -out extrainfo/$species2-$species1 -db extrainfo/$species1/$species1 -outfmt 6 -evalue 1e-5 -num_threads $pnum
fi

if [ ! -f "extrainfo/blastn_mapping_$species1-$species2.tab" ]; then
echo "summarize reciprocal blastN"
python __prepare_egio_blastn.py --species1 $species1 --blast1 extrainfo/$species1-$species2 --exonanno1 extrainfo/$species1.exon --species2 $species2 --blast2 extrainfo/$species2-$species1  --exonanno2 extrainfo/$species2.exon --orthogen $ortholog
rm extrainfo/$species1.exonfasta
rm extrainfo/$species2.exonfasta
rm extrainfo/$species1-$species2
rm extrainfo/$species2-$species1
fi

###=======================================
echo "detect orthologous exons and isoforms between species $species1 and $species2"
python __EGIO.py --orthog $ortholog --blastn extrainfo/blastn_mapping_$species1-$species2.tab --species1 $species1 --species2 $species2 --isocom1 extrainfo/$species1.tran --cdscom1 extrainfo/$species1.exon --isocom2 extrainfo/$species2.tran --cdscom2 extrainfo/$species2.exon --identhres $identity --coverthres $coverage --match $match --mismatch $mismatch --gap $gap --pnum $pnum

### now, the final results are stored in file: ExonGroup.txt and OrthoIso.txt
