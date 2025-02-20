#! /bin/bash

if [ $0 != "-bash" ]; then
    source $(dirname $0)/../filePaths.sh
else
    echo "Source the file 'filePaths.sh' in the github directory"
fi

condaEnvName=HintermannBoltHawkinsEtAl2025
# This line is to adapt the conda to the shell
source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
# Activate the conda environment
conda activate ${condaEnvName}

### ZEBRAFISH

# Get gtf:
wget https://zenodo.org/records/10283274/files/mergeOverlapGenesOfFilteredTranscriptsOfDanio_rerio.GRCz11.109_ExonsCDSOnly_UCSC.gtf.gz?download=1 -P ${genomeDirectory} -nc 

# Get the annotation colors 
wget https://raw.githubusercontent.com/lldelisle/scriptsForBoltEtAl2022/9644d5c4e3ea10b5dfae670fc08e25ddf2bbc78b/annotations/HoxD_Elements_mm10_colored.bed -P ${genomeDirectory} -nc

cd ${gitHubDirectory}/annotations

# subset in HoxD region and keep only protein coding
# As in the mouse, some isoforms of Hoxd3 are going over Hoxd4 (Hoxd3a-204 and Hoxd3a-201). There are also two lnpa isoforms.
# For simplification of figures we remove these isoforms.
gunzip -c ${genomeDirectory}/mergeOverlapGenesOfFilteredTranscriptsOfDanio_rerio.GRCz11.109_ExonsCDSOnly_UCSC.gtf.gz | awk -v OFS="\t" '$1=="chr9" && $4 < 2500000 && $5 > 1600000 {print $0}' | grep "protein_coding" | grep -v "hoxd3a-204\|hoxd3a-201\|lnpa-201"  > ${genomeDirectory}/danRer11_hoxDa_custom_protein_coding.gtf

# Get only one annotation per gene name and convert to bed
python ../scripts/fromgtfTobed12.py --output ${genomeDirectory}/danRer11_hoxDa_custom_protein_coding_details.bed --mergeTranscriptsAndOverlappingExons ${genomeDirectory}/danRer11_hoxDa_custom_protein_coding.gtf
# Make solid rectangles
awk -v OFS="\t" '{print $1,$2,$3,$4,"0",".",$2,$3,"163,163,163"}' ${genomeDirectory}/danRer11_hoxDa_custom_protein_coding_details.bed > ${genomeDirectory}/danRer11_hoxDa_custom_protein_coding.bed
# Get hox gene colors as copy of mouse colors
cat ${genomeDirectory}/HoxD_Elements_mm10_colored.bed | grep Hoxd | cut -f 4,9 > ${genomeDirectory}/Hoxd_colors.txt
echo -e "Hoxd3\t190,127,230" >> ${genomeDirectory}/Hoxd_colors.txt
awk -v OFS="\t" '/Hox/ {$1=tolower($1) "a"; print}' ${genomeDirectory}/Hoxd_colors.txt > ${genomeDirectory}/hoxDcolor_fish.txt

# Build final gene annotation files
# Get non-Hox genes
cat ${genomeDirectory}/danRer11_hoxDa_custom_protein_coding.bed | grep -v hox > danRer11_hoxDa_custom_protein_coding.bed
# Add colored Hox genes
join -1 4 -2 1 <(sort -k4,4 ${genomeDirectory}/danRer11_hoxDa_custom_protein_coding.bed) <(sort -k1,1 ${genomeDirectory}/hoxDcolor_fish.txt ) | awk -v OFS="\t" '{print $2,$3,$4,$1,$5,$6,$7,$8,$NF}' >> danRer11_hoxDa_custom_protein_coding.bed
# Make a version in which the hoxDa cluster is represented instead of the individual hoxda genes
hox_start=$(awk '$4 == "hoxd3a" {print $2}' danRer11_hoxDa_custom_protein_coding.bed)
hox_end=$(awk '$4 == "hoxd13a" {print $3}' danRer11_hoxDa_custom_protein_coding.bed)
awk -v OFS="\t" -v hox_start="$hox_start" -v hox_end="$hox_end" 'BEGIN {printed_hoxd = 0}{
    if (!($4 ~ /hoxd/)){
        $9="255,255,255";print;
    } else if (printed_hoxd == 0) {
        print "chr9", hox_start, hox_end, "hoxDa", "0", ".", hox_start, hox_end, "0,0,0";
        printed_hoxd = 1;
    }
}' danRer11_hoxDa_custom_protein_coding.bed > danRer11_hoxDa_custom_protein_coding_HoxDblack.bed
# Make a version of hoxda genes described to be expressed preaxially
awk -v OFS="\t" '
{
    if ($4 ~ /hoxd(3|4)a$/) {
        $9="163,163,163"
    }
    print
}' ${genomeDirectory}/danRer11_hoxDa_custom_protein_coding_details.bed > ${genomeDirectory}/danRer11_hoxDa_custom_protein_coding_details_tmp.bed 
awk -v OFS="\t" '
{
    if ($4 ~ /hoxd(9|10)a$/) {
        $9="150,51,201"
    }
    print
}' ${genomeDirectory}/danRer11_hoxDa_custom_protein_coding_details_tmp.bed > danRer11_hoxDa_custom_preax.bed
# Make a version of hoxda genes described to be expressed postaxially
awk -v OFS="\t" '{
    if ($4 ~ /^hoxd1([1-3])a$/) {
        $9="232,103,27"
    }
    print
}' ${genomeDirectory}/danRer11_hoxDa_custom_protein_coding_details_tmp.bed > danRer11_hoxDa_custom_postax.bed
# Generate label tracks
for f in danRer11_hoxDa_custom_protein_coding_HoxDblack.bed danRer11_hoxDa_custom_preax.bed; do 
    echo $f
    awk -v OFS="\t" '{print $1, $3, $3+1, $4}' $f > ${f/.bed/_start.bed}
done
