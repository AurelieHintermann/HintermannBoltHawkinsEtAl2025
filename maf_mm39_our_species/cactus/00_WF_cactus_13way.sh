#!/bin/bash

ssh cerebro
ml simr
sinteractive

ml cactus/2.6.7-gpu
ml pygenometracks

# https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md  
# https://github.com/ComparativeGenomicsToolkit/hal


### Config

# List of species to align.
# text file, col1 is common name (no space), col2 is assembly ID
# Aphyosemion_australe	GCA_006937985.1
# Danio_rerio	GCA_000002035.4

species_to_align="/n/projects/ah2717/vertebrate_alignments/vertebrate_alignments_genomes.txt"
path_to_genome_folder="/n/projects/ah2717/genomes/fasta/"

# the dirname of species_to_align is the working directory 
wdir=$(dirname $species_to_align)/

cactus_working="${wdir}working/"
cactus_jS="${wdir}jS/"
cactus_hal="${wdir}vertebrates_13way.hal"

# Dowload genomes
source /n/projects/ah2717/fish_alignments/functions.sh
cd $path_to_genome_folder

ln -s /n/analysis/indexes/mm39/mm39.fa Mus_musculus_GCA_000001635.9.fna
ln -s /n/analysis/indexes/danRer11/danRer11.fa Danio_rerio_GCA_000002035.4.fna
ln -s /n/analysis/indexes/hg38/hg38.fa Homo_sapiens_GCA_000001405.29.fna

while IFS=$'\t' read -r col1 col2; do
    # Define the expected file name
    genome_file="${col1}_${col2}.fna"
    
    # Check if the genome file already exists
    if [ ! -f "$genome_file" ]; then
        # Download the genome if it does not exist
        download_genome "$col2" "$col1"
        echo "Downloaded genome: $col2 | Saved as ${genome_file}"
    else
        # Notify that the genome file already exists
        echo "${genome_file} already exists. Not downloading."
    fi
done < $species_to_align

# Verify that each genome is soft masked (otherwise cactus takes too long)
while IFS=$'\t' read -r col1 col2; do
    genome_file="${col1}_${col2}.fna"
	if grep -q "[a-z]" "$genome_file"; then
    echo "Soft-masked: $genome_file"
  else
    echo "Not soft-masked: $genome_file"
  fi
done < $species_to_align

### Running cactus

cd $wdir

while IFS=$'\t' read -r col1 col2; do
    # Replace underscores with spaces
    sp_no_underscore=$(echo "$col1" | sed 's/_/ /g')
    # Write to the output file
    echo "$sp_no_underscore" >> "${species_to_align/.txt/_timetreeIn.txt}"
done < $species_to_align

# Make a tree that you call ${species_to_align/.txt/_timetreeIn_noBrL.nwk}
# You can upload "${species_to_align/.txt/_timetreeIn.txt}" in https://timetree.org/, then download the nwk file and remove branch lengths
# (time branch lengths don't work with cactus)

sbatch /n/projects/ah2717/fish_alignments/scripts/02_cactus.sh $species_to_align $path_to_genome_folder $cactus_working $cactus_jS $cactus_hal

######
### hal tools (post alignment)
######

# https://github.com/ComparativeGenomicsToolkit/hal
# https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md


halValidate vertebrates_13way.hal
# Returned: File valid
halStats vertebrates_13way.hal

# $1 = jobStore folder
# $2 = halFile
# $3 = outputMAF
# $4 = --refGenome
# $5 = --chunkSize

sbatch /n/projects/ah2717/fish_alignments/scripts/05_run_hal2maf.sh $cactus_jS $cactus_hal ${cactus_hal/.hal/_projected_Dr.maf} Danio_rerio 500000
sbatch /n/projects/ah2717/fish_alignments/scripts/05_run_hal2maf.sh ${cactus_jS/jS/jS2} $cactus_hal ${cactus_hal/.hal/_projected_Mm.maf} Mus_musculus 500000

# Project only chr of interest
cactus-hal2maf $cactus_jS $cactus_hal ${cactus_hal/.hal/_projected_Dr_chr9.maf} --refGenome Danio_rerio --chunkSize 500000 --refSequence chr9 --noAncestors
cactus-hal2maf ${cactus_jS/jS/jS2} $cactus_hal ${cactus_hal/.hal/_projected_Mm_chr2.maf} --refGenome Mus_musculus --chunkSize 500000 --refSequence chr2 --noAncestors

