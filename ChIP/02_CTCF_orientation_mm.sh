#!/bin/bash

if [ $0 != "-bash" ]; then
    source $(dirname $0)/../filePaths.sh
else
    echo "Source the file 'filePaths.sh' in the github directory"
fi

pathToFasta="${genomeDirectory}/fasta/"
ChIP_CTCF_mm="${GEODirectory}/ChIP/CTCF_ChIP_mm_E105_PT_rep1.narrowPeak.gz"
# Generate CTCF motif tables
generateNewMotifPrediction_forMouse=yes

condaEnvName=HintermannBoltEtAl2023
# This line is to adapt the conda to the shell
source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
# Activate the conda environment
conda activate ${condaEnvName}

echo "This script is run in $PWD. It will generate a lot of intermediate files which can be removed afterward."

source "${gitHubDirectory}/scripts/functions.sh"

# Make output directory
mkdir -p "${gitHubDirectory}/ChIP/CTCF_orientation"

### Mouse CTCF peaks
short_name=$(basename ${ChIP_CTCF_mm/.narrowPeak.gz/})

if [ "$generateNewMotifPrediction_forMouse" = "yes" ]; then

    # Get the chr2 fasta:
    wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/chromosomes/chr2.fa.gz -P ${pathToFasta} -nc
    if [ ! -e ${pathToFasta}/chr2.fa ]; then
        gunzip -k ${pathToFasta}/chr2.fa.gz
    fi
    # Usage: <peaks.narrowPeak.gz, chr, start, end, genome.fasta>
    # Description: this function takes the narrowPeak as input, filters peaks in a region of interest, extends each peak by 100 bp each side and extract the sequence of each peak. 
    extendPeaksAndGetFasta "$ChIP_CTCF_mm" chr2 73700000 75800320 ${pathToFasta}/chr2.fa

    printf "Upload %s_extended.fa on the website of http://insulatordb.uthsc.edu in CTCFBS Prediction Tool.\n" $PWD/$short_name
    read -p "Save table output (without header) as ${short_name}_extended_insdb_output.txt and press any key to resume."

fi 

# Usage: reformat_insulatordb_table <*_extended_insdb_output.txt> <*_noEMBL.bed>
# Description: this function takes the insdb table as input, removes the EMBL motives (because very short), put . if motif score is negative and keeps the best motif for each peak.
# Output: <*_noEMBL.bed>

reformat_insulatordb_table ${short_name}_extended_insdb_output.txt ${short_name}_extended_noEMBL.bed

# Function to create a bed with rgb field corresponding to motif orientation of CTCF
# Usage: makeRgbField <*extended_noEMBL.bed> [inverted]

makeRgbField ${short_name}_extended_noEMBL.bed

# Make a table where each motif is associated to the peak
bedtools intersect -a ${short_name}_extended_noEMBL_colored.bed -b ${ChIP_CTCF_mm} -wa -wb > ${short_name}_extended_noEMBL_colored_withPeakInfo.bed

# Select only peaks above a peak score (here 500)
peakScoreThreshold=500
awk -F "\t" -v OFS="\t" -v t=$peakScoreThreshold '
{
    if ($14 > t){
        # Print only the bed9 with motif
        print $1, $2, $3, $4, $5, $6, $7, $8, $9
    }
}' ${short_name}_extended_noEMBL_colored_withPeakInfo.bed > "${gitHubDirectory}/ChIP/CTCF_orientation/${short_name}_extended_noEMBL_colored_peakScorePlus500.bed"
