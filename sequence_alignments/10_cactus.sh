#!/bin/bash

# cactus version 2.6.7
# samtools version 1.20

######
### 13-way alignment
######

path_to_genome_folder="/n/projects/ah2717/genomes/fasta/"
cd $path_to_genome_folder

species_to_align="${gitHubDirectory}/sequence_alignments/cactus/13way/vertebrate_alignments_genomes.txt"

source "${gitHubDirectory}/sequence_alignments/11_download_genomes.sh"

### Dowload genomes

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

### Get a species tree

# species tree
# Save species list to draw nice trees in timeTree
# http://timetree.org/
	
while IFS=$'\t' read -r col1 col2; do
    # Replace underscores with spaces
    sp_no_underscore=$(echo "$col1" | sed 's/_/ /g')
    # Write to the output file
    echo "$sp_no_underscore" >> "${species_to_align/.txt/_timetreeIn.txt}"
done < $species_to_align

# Build a Timetree >  Load a List of Species
# Leucoraja erinacea (replaced with Amblyraja radiata)
# "to Newick File", saved as vertebrate_alignments_genomes_timetreeIn.nwk
# delete branch lenght (time tree has time branch lengths instead of probability branch lengths and it makes cactus crash)
# final file is the input for cactus:
	# vertebrate_alignments_genomes_timetreeIn_noBrL.nwk

### Run cactus

# the dirname of species_to_align is the working directory 
wdir=$(dirname $species_to_align)/
cd $wdir

cactus_working="${wdir}working/"
cactus_jS="${wdir}jS/"
cactus_hal="${wdir}vertebrates_13way.hal"

sbatch "${gitHubDirectory}/sequence_alignments/12_cactus.sh" $species_to_align $path_to_genome_folder $cactus_working $cactus_jS $cactus_hal
# The .hal file is available on Zenodo at https://zenodo.org/records/14825508

# Verify hal
halValidate $cactus_hal

# Project on chr of interest
cactus-hal2maf ${cactus_jS/jS/jS2} $cactus_hal ${cactus_hal/.hal/_projected_Dr_chr9.maf} --refGenome Danio_rerio --chunkSize 500000 --refSequence chr9 --noAncestors
# Cannot have twice the same "jS" file
cactus-hal2maf ${cactus_jS/jS/jS3} $cactus_hal ${cactus_hal/.hal/_projected_Mm_chr2.maf} --refGenome Mus_musculus --chunkSize 500000 --refSequence chr2 --noAncestors

######
### 3-way alignment on HoxD and 5DOM
######

### Run cactus

# Use the whole mouse chr2 so that we can use mouse annotations

# extract HoxD and 5DOM region in zebrafish and gar genomes to ensure we are looking at homologous regions when projecting on mouse chr2:73640000-74630000
samtools faidx ${path_to_genome_folder}Danio_rerio_GCA_000002035.4.fna chr9:1900000-2400000 > ${path_to_genome_folder}Danio_rerio_GCA_000002035.4_HoxD5DOM.fna
samtools faidx ${path_to_genome_folder}Lepisosteus_oculatus_GCF_000242695.1.fna NC_023190.1:15506000-15815550 > ${path_to_genome_folder}Lepisosteus_oculatus_GCF_000242695.1_HoxD5DOM.fna

# 
species_to_align="${gitHubDirectory}/sequence_alignments/cactus/3way_HoxD_5DOM/mm_zf_lo_alignments_genomes.txt"

# the dirname of species_to_align is the working directory 
wdir=$(dirname $species_to_align)/
cd $wdir

cactus_working="${wdir}working/"
cactus_jS="${wdir}jS/"
cactus_hal="${wdir}vertebrates_13way.hal"

sbatch "${gitHubDirectory}/sequence_alignments/12_cactus.sh" $species_to_align $path_to_genome_folder $cactus_working $cactus_jS $cactus_hal
# The .hal file is available on Zenodo at https://zenodo.org/records/14825508

# Verify hal
halValidate $cactus_hal

# Project on chr of interest
cactus-hal2maf ${cactus_jS/jS/jS2} $cactus_hal ${cactus_hal/.hal/_projected_Mm.maf} --refGenome Mus_musculus --chunkSize 500000 --noAncestors


