#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 10G # The memory needed depends on the size of the genome and the size of fastqs
#SBATCH --cpus-per-task 1 # This allows to speed the mapping part of HiCUP
#SBATCH --time 48:00:00 # This depends on the size of the fastqs
#SBATCH --array=2-9 # Put here the rows from the table that need to be processed in the table
#SBATCH --job-name chaining_maffing # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /scratch/ldelisle/maf/ # This directory must exist, this is where will be the error and out files
## Specific to jed:
#SBATCH --qos serial

scriptDirectory=$1
filePathForTable=$2
fastaS=$3

module purge
export PATH=$scriptDirectory:$PATH

# The parameters for chainFar were taken from
# https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsGalGal6/
# and the header of
# https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsGalGal6/mm39.galGal6.all.chain.gz
chainFar="-minScore=5000 -linearGap=loose"

# Get the assembly from the table
assemblyQ=$(cat ${filePathForTable} | gawk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}')
# Get chromosome name and assembly name
chromNameS=$(basename $fastaS .fa)
assemblyS=$(basename $(dirname $fastaS))

globalname=${assemblyS}.${chromNameS}.${assemblyQ}

# In the UCSC pipeline,
# The code for this step is
# https://github.com/ENCODE-DCC/kentUtils/blob/master/src/hg/utils/automation/doBlastzChainNet.pl

# To help me to find the commands
# I also used
# https://genomewiki.ucsc.edu/images/9/93/RunLastzChain_sh.txt

mkdir -p chain
for psl_file in psl/${globalname}.*.psl.gz; do
    echo $psl_file
    name=$(basename $psl_file .psl.gz)
    chromNameQ=$(echo $name | awk -v as=$assemblyQ '{split($0, a, as"."); print a[2]}')
    echo $assemblyQ $chromNameQ
    zcat $psl_file | axtChain -psl -faQ -faT ${chainFar} stdin ${fastaS} ${assemblyQ}/${chromNameQ}.fa chain/${name}.chain
done

mkdir -p temp/
chainMergeSort chain/${globalname}.*.chain > temp/${globalname}.chain
chainPreNet temp/${globalname}.chain  ${assemblyS}/${assemblyS}.sizes ${assemblyQ}/${assemblyQ}.sizes temp/${globalname}.pre.chain
chainNet temp/${globalname}.pre.chain -minSpace=1 ${assemblyS}/${assemblyS}.sizes ${assemblyQ}/${assemblyQ}.sizes stdout /dev/null | netSyntenic stdin temp/${globalname}.net
# On https://genomewiki.ucsc.edu/images/9/93/RunLastzChain_sh.txt
# I should do:
# netClass -noAr temp/${globalname}.net ${assemblyS} ${assemblyQ} temp/${globalname}.netClass.net
# netToAxt temp/${globalname}.netClass.net temp/${globalname}.pre.chain ${assemblyS}/${assemblyS}.2bit ${assemblyQ}/${assemblyQ}.2bit stdout | axtSort stdin temp/${globalname}.axt

# On the doBlastzChainNet.pl there is a netFilter step
# https://github.com/ENCODE-DCC/kentUtils/blob/3e7ba71fae2070be32d59bd5b5ba70a6c9aa9c06/src/hg/utils/automation/doBlastzChainNet.pl#L1559


netToAxt temp/${globalname}.net temp/${globalname}.pre.chain ${assemblyS}/${assemblyS}.2bit ${assemblyQ}/${assemblyQ}.2bit stdout | axtSort stdin temp/${globalname}.axt

mkdir -p maf_per_species
axtToMaf temp/${globalname}.axt ${assemblyS}/${assemblyS}.sizes ${assemblyQ}/${assemblyQ}.sizes maf_per_species/${globalname}.maf -tPrefix=${assemblyS}. -qPrefix=${assemblyQ}.
