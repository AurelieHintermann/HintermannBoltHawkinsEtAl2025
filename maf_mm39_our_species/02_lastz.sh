#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 10G # The memory needed depends on the size of the genome and the size of fastqs
#SBATCH --cpus-per-task 1 # This allows to speed the mapping part of HiCUP
#SBATCH --time 48:00:00 # This depends on the size of the fastqs
#SBATCH --array=2-9 # Put here the rows from the table that need to be processed in the table
#SBATCH --job-name lastz # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /scratch/ldelisle/maf/ # This directory must exist, this is where will be the error and out files
## Specific to jed:
#SBATCH --qos serial

scriptDirectory=$1
filePathForTable=$2
fastaS=$3

module purge
export PATH=$scriptDirectory:$PATH

# The parameters for lastzFar were taken from
# https://genomewiki.ucsc.edu/index.php/Mm39_35-way_conservation_lastz_parameters
# As well as doc from 
# https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsGalGal6/
# and the header of
# https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsGalGal6/mm39.galGal6.all.chain.gz

# lastzFar="E=30 H=2000 K=3000 L=3000 M=0 O=400 Q=${scriptDirectory}/../default.q T=1 Y=9400"
lastzFar="H=2000"

# Get the assembly from the table
assemblyQ=$(cat ${filePathForTable} | gawk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}')
# Get chromosome name and assembly name
chromNameS=$(basename $fastaS .fa)
assemblyS=$(basename $(dirname $fastaS))

fastaQs=$(find $assemblyQ -name "*.fa")

# In the UCSC pipeline,
# The code for this step is
# https://github.com/ENCODE-DCC/kentUtils/blob/master/src/hg/utils/automation/blastz-run-ucsc
# But there is a partition step to cut the chrom seq into chunks
# The master script is 
# https://github.com/ENCODE-DCC/kentUtils/blob/master/src/hg/utils/automation/doBlastzChainNet.pl#L484

# To help me to find the commands
# I also used
# https://genomewiki.ucsc.edu/images/9/93/RunLastzChain_sh.txt

mkdir -p raw
mkdir -p psl

for fastaQ in $fastaQs; do
    chromNameQ=$(basename $fastaQ .fa)
    name=${assemblyS}.${chromNameS}.${assemblyQ}.${chromNameQ}
    echo $name

    if [ ! -e raw/${name}.lav ]; then
        lastz-1.04.00 $fastaS $fastaQ $lastzFar > raw/${name}.lav
    fi
    if [ ! -e psl/${name}.psl.gz ]; then
        lavToPsl raw/${name}.lav stdout \
            | gzip -c > psl/${name}.psl.gz
    fi
done
