#!/bin/bash

# The goal of this bash file is to generate annotations/danRer11_hoxDa_elements.bed
# Which contains annotations/mm39_HoxD_elements.bed which are homologous in danRer11
# And to generate bed files with conserved elements
# Between mouse and zebrafish

if [ $0 != "-bash" ]; then
    source $(dirname $0)/../filePaths.sh
else
    echo "Source the file 'filePaths.sh' in the github directory"
fi

pathToFasta="${genomeDirectory}/fasta/"

cd $gitHubDirectory

condaEnvName=HintermannBoltEtAl2023
# This line is to adapt the conda to the shell
source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
# Activate the conda environment
conda activate ${condaEnvName}

# Create output directory
mkdir -p sequence_alignments/conserved_elements

# Define regions to align
mm39_HoxD_fasta=("chr2" 73600000 75670000)
danRer11_HoxDa_fasta=("chr9" 1639000 2400000)

# Extract fasta for regions of interest
# mm39
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz -P ${pathToFasta} -nc
if [ ! -e ${pathToFasta}/mm39.fa ]; then
    gunzip -k ${pathToFasta}/mm39.fa.gz
fi
echo -e "${mm39_HoxD_fasta[0]}\t${mm39_HoxD_fasta[1]}\t${mm39_HoxD_fasta[2]}" > ${pathToFasta}/mm39_HoxD.bed
bedtools getfasta -fi ${pathToFasta}/mm39.fa  -bed ${pathToFasta}/mm39_HoxD.bed > ${pathToFasta}/mm39_HoxD.fa

# danRer11
if [ ! -e ${pathToFasta}/danRer11_noalt_noUn.fa ]; then
    bash prepare/01_get_fasta_and_gtf.sh
fi
echo -e "${danRer11_HoxDa_fasta[0]}\t${danRer11_HoxDa_fasta[1]}\t${danRer11_HoxDa_fasta[2]}" > ${pathToFasta}/danRer11_HoxDa.bed
bedtools getfasta -fi ${pathToFasta}/danRer11_noalt_noUn.fa -bed ${pathToFasta}/danRer11_HoxDa.bed > ${pathToFasta}/danRer11_HoxDa.fa

# Align fish to mouse
blastn -subject ${pathToFasta}/mm39_HoxD.fa -query ${pathToFasta}/danRer11_HoxDa.fa -task dc-megablast -outfmt 6 > ${pathToFasta}/alignment_mm_zf_discontinuousMegablast.txt
# Align mouse to fish 
blastn -subject ${pathToFasta}/danRer11_HoxDa.fa -query ${pathToFasta}/mm39_HoxD.fa -task dc-megablast -outfmt 6 > ${pathToFasta}/alignment_zf_mm_discontinuousMegablast.txt

# Add header
for f in ${pathToFasta}/*discontinuousMegablast.txt; do 
    echo -e "qaccver\tsaccver\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" \
    | cat - $f > ${pathToFasta}/tmp.txt
    mv ${pathToFasta}/tmp.txt $f
done

for f in ${pathToFasta}/*discontinuousMegablast.txt; do 
    awk -F "\t" -v OFS="\t" '
NR > 1{
    print substr($1,1,index($1,":")-1), $7, $8, substr($2,1,index($2,":")-1), $9, $10
}' $f > ${f/.txt/_coordinate_correspondance.txt}
done

# Keep only sequences with reciprocal hits and sort by mouse coordinates
awk -F "\t" -v OFS="\t" -v mm_start="${mm39_HoxD_fasta[1]}" -v zf_start="${danRer11_HoxDa_fasta[1]}" '
NR==FNR{
    a[$4,$6,$5,$1,$3,$2]
    next
}
NR>1 && ($1,$2,$3,$4,$5,$6) in a {
    print $1,zf_start+$2,zf_start+$3,$4,$5+mm_start,$6+mm_start
}' ${pathToFasta}/alignment_zf_mm_discontinuousMegablast_coordinate_correspondance.txt ${pathToFasta}/alignment_mm_zf_discontinuousMegablast_coordinate_correspondance.txt | sort -k5,5n > ${pathToFasta}/discontinuousMegablast_reciprocal_hits_tmp.txt
# Add header and attribute unique name to each conserved sequence
awk -F "\t" -v OFS="\t" '
NR==1{
    print "zf_chr", "zf_start", "zf_end", "mm_chr", "mm_start", "mm_end", "CS_ID"
    next
}
NR==FNR{
    print $0,"CS_"FNR-1
}' ${pathToFasta}/discontinuousMegablast_reciprocal_hits_tmp.txt > sequence_alignments/conserved_elements/discontinuousMegablast_reciprocal_hits.txt
rm ${pathToFasta}/discontinuousMegablast_reciprocal_hits_tmp.txt

# Write conserved elements bed file for each species
awk -F "\t" -v OFS="\t" '
NR == 1{
    print "#chr", "start", "end", "CS_ID"
}
NR > 1{
    print $4,$6,$5,$7
}' sequence_alignments/conserved_elements/discontinuousMegablast_reciprocal_hits.txt > sequence_alignments/conserved_elements/mm_conserved_elements_mm_zf.bed
awk -F "\t" -v OFS="\t" '
NR == 1{
    print "#chr", "start", "end", "CS_ID"
}
NR > 1{
    print $1,$2,$3,$7
}' sequence_alignments/conserved_elements/discontinuousMegablast_reciprocal_hits.txt > sequence_alignments/conserved_elements/zf_conserved_elements_mm_zf.bed

## Annotate fish CS which were described as enhancers in mouse
# Find mouse CS which overlap mouse elements
bedtools intersect -a annotations/mm39_HoxD_elements.bed -b sequence_alignments/conserved_elements/mm_conserved_elements_mm_zf.bed -wa -wb | awk -F "\t" -v OFS="\t" '{print $4,$13,$9}' > ${pathToFasta}/tmp1.bed
# Find correspond fish Cs
join -1 4 -2 2 -o '1.1 1.2 1.3 1.4 2.1 2.3' -t $'\t' <(sort -k4 sequence_alignments/conserved_elements/zf_conserved_elements_mm_zf.bed) <(sort -k2 ${pathToFasta}/tmp1.bed) > ${pathToFasta}/tmp2.bed
# Merge regions by element names
awk -F "\t" -v OFS="\t" '
!($5 in min || $5 in max) ($2 < min[$5]) {
    min[$5] = $2
}
($3 > max[$5]) {
    max[$5] = $3
}
END {
    for (i in min){
        print $1, min[i], max[i], i, "0", ".",min[i], max[i],$6
    }
}' ${pathToFasta}/tmp2.bed > annotations/danRer11_hoxDa_elements.bed
