#!/bin/bash

#SBATCH -o err_out/slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e err_out/slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --partition gpu # Tell SLURM to look for nodes in the 'GPU' partition
#SBATCH --gres=gpu:a100:1 # Request 1 GPU
#SBATCH -c 40
#SBATCH --mem=256G
#SBATCH --job-name=cactus
#SBATCH --time=6-00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ahintermann@stowers.org

ml cactus/2.4.0-gpu

export OMP_NUM_THREADS=1

# Meaning of variables
# $1=species_to_align
# $2=path_to_genome_folder
# $3=cactus_working
# $4=cactus_jS (job store)
# $5=cactus_hal (output)

# Record the start time
start_time=$(date +%s)

cd $(dirname $1)

### Format cactus input file (= seq file)
# <seqFile> : Text file containing the locations of the input sequences as well as their phylogenetic tree. 

# write tree as line 1
cat ${1/.txt/_timetreeIn_noBrL.nwk} > ${1/.txt/_tree.txt} 

# write species names (col1) and path to genome (col2)
# This file is needed as cactus input
while IFS=$'\t' read -r col1 col2; do
    # Define the expected file name
    genome_file="${col1}_${col2}.fna"
    
    # Check if the genome file already exists
    if [ ! -f "$2$genome_file" ]; then
        # Download the genome if it does not exist
        echo "${genome_file} does not exist."
    else
        # Notify that the genome file already exists
        echo "${genome_file} exists."
        echo -e "${col1}\t$2${genome_file}" >> ${1/.txt/_tree.txt} 
    fi
done < $1

### run cactus
# make working dir
mkdir -p $3

cactus --gpu 1  \
    --clean always \
    --workDir $3 \
    $4 \
    ${1/.txt/_tree.txt} \
    $5
#    --restart \

### write commands in .out file
echo "

### Format cactus input file (= seq file)
# <seqFile> : Text file containing the locations of the input sequences as well as their phylogenetic tree. 

# write tree as line 1
cat ${1/.txt/_timetreeIn_noBrL.nwk} > ${1/.txt/_tree.txt} 

# write species names (col1) and path to genome (col2)
# This file is needed as cactus input
while IFS=$'\t' read -r col1 col2; do
    # Define the expected file name
    genome_file="${col1}_${col2}.fna"
    
    # Check if the genome file already exists
    if [ ! -f "$2$genome_file" ]; then
        # Download the genome if it does not exist
        echo "${genome_file} does not exist."
    else
        # Notify that the genome file already exists
        echo "${genome_file} exists."
        echo -e "${col1}\t$2${genome_file}" >> ${1/.txt/_tree.txt} 
    fi
done < $1

### run cactus
# make working dir
mkdir -p $3

cactus --gpu 1  \
    --clean always \
    --workDir $3 \
    $4 \
    ${1/.txt/_tree.txt} \
    $5
#    --restart \

"
# Record the end time
end_time=$(date +%s)
# Calculate the elapsed time
elapsed_time=$(( (end_time - start_time) / 60 ))
echo "Job completed in ${elapsed_time} minutes."