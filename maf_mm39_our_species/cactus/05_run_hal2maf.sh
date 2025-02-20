#!/bin/bash

#SBATCH --job-name=hal2maf
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH -o err_out/slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e err_out/slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ahintermann@stowers.org

ml cactus/2.4.0-gpu

# Record the start time
start_time=$(date +%s)

# $1 = jobStore folder
# $2 = halFile
# $3 = outputMAF
# $4 = --refGenome
# $5 = --chunkSize

cactus-hal2maf $1 $2 $3 --refGenome $4 --chunkSize $5

# Record the end time
end_time=$(date +%s)
# Calculate the elapsed time
elapsed_time=$(( (end_time - start_time) / 60 ))
echo "Job completed in ${elapsed_time} minutes."