#!/bin/bash

#SBATCH -o slurm-%x-%A.out # Template for the std output of the job uses the job name and the job id
#SBATCH -e slurm-%x-%A.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 50G # The memory needed depends on the size of the genome and the size of fastqs
#SBATCH --cpus-per-task 1 # This only uses a single CPU
#SBATCH --time 12:00:00 # This depends on the number of cool files
#SBATCH --job-name TAD_calling # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /scratch/ldelisle/HintermannBolt/TAD_calling/ # This directory must exist, this is where will be the error and out files
#SBATCH --qos=serial

gitDir=$1
GEODirectory=$2

# This script call TADs

filePathForTADcallingOutput="${gitDir}/HiCandcHiC/TADcalling/"
condaEnvName=HintermannBoltEtAl2023
# This line is to adapt the conda to the shell
source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
# Activate the conda environment
conda activate ${condaEnvName}

mkdir -p "${filePathForTADcallingOutput}"

# call TADs on mouse UGS Capture-HiC-seq data
hicFindTADs --matrix ${GEODirectory}/cHiC/CHiC_mm_E185_MUGS.10kb.cool --outPrefix CHiC_mm_E185_MUGS.10kb_240kb --correctForMultipleTesting fdr --minDepth 240000 --maxDepth 480000 --step 480000 --minBoundaryDistance 200000

# Put the domains in the output directory
cp CHiC_mm_E185_MUGS.10kb_240kb_domains.bed ${filePathForTADcallingOutput}/

# call TADs on zebrafish whole embryo HiC-seq data
hicFindTADs --matrix ${GEODirectory}/HiC/HiC_hpf2448merge_WE_10kb.cool --outPrefix HiC_hpf2448merge_WE_10kb_70kb --chromosomes "chr9" --correctForMultipleTesting fdr --minDepth 70000 --maxDepth 140000 --step 140000 --minBoundaryDistance 60000

# Put the domains in the output directory
cp HiC_hpf2448merge_WE_10kb_70kb_domains.bed ${filePathForTADcallingOutput}/
