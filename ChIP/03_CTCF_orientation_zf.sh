#!/bin/bash

if [ $0 != "-bash" ]; then
    source $(dirname $0)/../filePaths.sh
else
    echo "Source the file 'filePaths.sh' in the github directory"
fi

pathToFasta="${genomeDirectory}/fasta/"
ChIP_CTCF_zf_prefix="${GEODirectory}/ChIP/Franke_CTCF_ChIP_"
short_name=$(basename ${ChIP_CTCF_zf_prefix})
# Generate CTCF motif tables
generateNewMotifPrediction_forZf=yes

condaEnvName=HintermannBoltEtAl2023
# This line is to adapt the conda to the shell
source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
# Activate the conda environment
conda activate ${condaEnvName}

echo "This script is run in $PWD. It will generate a lot of intermediate files which can be removed afterward."

source "${gitHubDirectory}/scripts/functions.sh"

# Make output directory
mkdir -p "${gitHubDirectory}/ChIP/CTCF_orientation"


if [ "$generateNewMotifPrediction_forZf" = "yes" ]; then

    # Get the danRer11 fasta:
    if [ ! -e ${pathToFasta}danRer11_noalt_noUn.fa ]; then
        bash prepare/01_get_fasta_and_gtf.sh
    fi
    
    # Get fasta file for extended CTCF peaks in region of interest
    for rep in ${ChIP_CTCF_zf_prefix}*.narrowPeak.gz; do
        # Usage: <peaks.narrowPeak.gz, chr, start, end, genome.fasta>
        # Description: this function takes the narrowPeak as input, filters peaks in a region of interest, extends each peak by 100 bp each side and extract the sequence of each peak. 
        extendPeaksAndGetFasta "$rep" chr9 1600000 2500000 ${pathToFasta}/danRer11_noalt_noUn.fa
    done

    printf "Upload each %s*extended.fa on the website of http://insulatordb.uthsc.edu in CTCFBS Prediction Tool.\n" $PWD/${short_name}
    read -p "Save each output table (without header) as extended_insdb_output.txt and press any key to resume."
fi

for tbl in ${short_name}*extended_insdb_output.txt; do
    # Usage: reformat_insulatordb_table <*_extended_insdb_output.txt> <*_noEMBL.bed>
    # Description: this function takes the insdb table as input, removes the EMBL motives (because very short), put . if motif score is negative and keeps the best motif for each peak.
    # Output: <*_noEMBL.bed>
    reformat_insulatordb_table "$tbl" "${tbl/_insdb_output.txt/_noEMBL.bed}"
done

# Keep peaks if they are present in the two replicates of at least one stage

# Get the peaks which are in both replicates, for each stage
intersect_peaks=""
for time in 24 48; do
    # Concatenate peaks of stage replicates in one file
    cat ${short_name}${time}hpf_rep[12]_extended_noEMBL.bed | sort -k1,1  -k2,2n > ${short_name}${time}hpf_allReps_allPeaks_noEMBL.bed
    # Description: this function takes the concatenate file of all replicates and the number of replicates as input. It keeps only motives found in all replicates, and keeps the one with the best score.
    # Usage: keepMotivesRepIntersection <cat_sorted_allPeaks_allReps.bed> <nbReps> <cat_sorted_allReps.bed>
    keepMotivesRepIntersection ${short_name}${time}hpf_allReps_allPeaks_noEMBL.bed 2 ${short_name}${time}hpf_allReps_noEMBL.bed
    intersect_peaks="$intersect_peaks ${short_name}${time}hpf_allReps_noEMBL.bed"
done

# Concatenate peaks of stage intersection
cat $intersect_peaks | sort -k1,1  -k2,2n -k3,3 > ${short_name}2448hpf_allReps_noEMBL.bed

# Get peaks which are in at least one stage
# Description: this function takes the concatenate file of all replicates and the number of replicates as input. For each motif, it keeps the one with the best score.
# Usage: keepMotivesRepUnion <cat_sorted_allPeaks_allReps.bed> <cat_sorted_atLeastOneRep.bed>
keepMotivesRepUnion ${short_name}2448hpf_allReps_noEMBL.bed ${short_name}allReps_noEMBL.bed

# Function to create a bed with rgb field corresponding to motif orientation of CTCF
# Usage: makeRgbField <*extended_noEMBL.bed> [inverted]
makeRgbField ${short_name}allReps_noEMBL.bed inverted
cp ${short_name}allReps_noEMBL_colored.bed "${gitHubDirectory}/ChIP/CTCF_orientation/"
