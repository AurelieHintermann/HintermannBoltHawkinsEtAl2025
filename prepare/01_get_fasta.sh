#! /bin/bash

if [ $0 != "-bash" ]; then
    source $(dirname $0)/../filePaths.sh
else
    echo "Source the file 'filePaths.sh' in the github directory"
fi

condaEnvName=HintermannBoltHawkinsEtAl2025
# This line is to adapt the conda to the shell
source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
# Activate the conda environment
conda activate ${condaEnvName}

mkdir -p ${genomeDirectory}

# First get fasta:
cd ${genomeDirectory}
mkdir fasta
cd fasta

# Get mm39 genome fasta
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz -nc
if [ ! -e mm39.fa ]; then
    gunzip -k mm39.fa.gz
fi

# Make danRer11_noalt
# Get the danRer11 genome fasta
wget https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/danRer11.fa.gz -nc
if [ ! -e danRer11.fa ]; then
    gunzip -k danRer11.fa.gz
fi
# Get the chrom size:
wget https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/danRer11.chrom.sizes -nc 
# Get the list of chromosome to keep
cut -f 1 danRer11.chrom.sizes | grep -v "_alt$" > noalt.txt
# Generate the noalt fasta
if [ ! -e danRer11_noalt.fa ]; then
    seqtk subseq -l 50 danRer11.fa noalt.txt > danRer11_noalt.fa
fi

# For the homology sequence analysis and for CTCF motif, use no alt and no chrUn:
cut -f 1 danRer11.chrom.sizes | grep -v "_alt$" | grep -v "^chrUn" > noalt_noUn.txt
if [ ! -e danRer11_noalt_noUn.fa ]; then
    seqtk subseq -l 50 danRer11.fa noalt_noUn.txt > danRer11_noalt_noUn.fa
fi

samtools faidx danRer11_noalt_noUn.fa
