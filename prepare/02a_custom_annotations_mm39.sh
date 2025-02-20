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

cd $gitHubDirectory

## Get files and scripts
# Get mm39 gtf
wget "https://zenodo.org/record/7510797/files/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm39.108_ExonsCDSOnly_UCSC.gtf.gz" -P ${genomeDirectory} -nc 
# Get the mm39 chr2 fasta
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/chromosomes/chr2.fa.gz -P ${genomeDirectory} -nc
# Get the annotation colors from Bolt et al 2022
wget https://raw.githubusercontent.com/lldelisle/scriptsForBoltEtAl2022/9644d5c4e3ea10b5dfae670fc08e25ddf2bbc78b/annotations/HoxD_Elements_mm10_colored.bed -P ${genomeDirectory} -nc
# Get script to convert gtf to bed from Guerrero et al. 2020
wget "https://raw.githubusercontent.com/lldelisle/scriptsForFernandezGuerreroEtAl2020/e399194b40999f8fd24c4370561813844039150a/RNAseq/fromgtfTobed12.py" -P scripts -nc
# Get the script for in silico PCR from @egonozer
wget https://raw.githubusercontent.com/egonozer/in_silico_pcr/7dd1ee2d3361d42bf4111b8bc77d8120e6d097f7/in_silico_PCR.pl -P scripts -nc

## Get a simplified gene annotation file
cd annotations
gunzip -k ${genomeDirectory}/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm39.108_ExonsCDSOnly_UCSC.gtf.gz
# Keep only protein coding transcripts, and keep one annotation per Hoxd gene
grep protein_coding ${genomeDirectory}/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm39.108_ExonsCDSOnly_UCSC.gtf | grep -v "Hoxa3-201" | grep -v "Hoxd3-203" > ${genomeDirectory}/mm39_custom_protein_coding.gtf
# Retain only entries within our region of interest
cat ${genomeDirectory}/mm39_custom_protein_coding.gtf | awk '$1=="chr2" && $4 < 77000000 && $5 > 71500000{print}' > ${genomeDirectory}/mm39_HoxD_custom_protein_coding.gtf
# Get only one annotation per gene name and convert to bed
python ${gitHubDirectory}/scripts/fromgtfTobed12.py --output mm39_HoxD_custom_protein_coding_details.bed --mergeTranscriptsAndOverlappingExons ${genomeDirectory}/mm39_HoxD_custom_protein_coding.gtf
# Convert bed12 to bed9 with color
awk -v OFS="\t" '{print $1,$2,$3,$4,"0",".",$2,$3,"163,163,163"}' mm39_HoxD_custom_protein_coding_details.bed > ${genomeDirectory}/mm39_HoxD_custom_protein_coding_tmp.bed
# Get colors for Hoxd genes
cat ${genomeDirectory}/HoxD_Elements_mm10_colored.bed | grep Hoxd | cut -f 4,9 > ${genomeDirectory}/Hoxd_colors.txt
echo -e "Hoxd3\t190,127,230" >> ${genomeDirectory}/Hoxd_colors.txt

# Build final gene annotation files
# Get non-Hox genes
cat ${genomeDirectory}/mm39_HoxD_custom_protein_coding_tmp.bed | grep -v Hox > mm39_HoxD_custom_protein_coding.bed
# Add colored Hox genes
join -1 4 -2 1 <(sort -k4,4 ${genomeDirectory}/mm39_HoxD_custom_protein_coding_tmp.bed) <(sort -k1,1 ${genomeDirectory}/Hoxd_colors.txt ) | awk -v OFS="\t" '{print $2,$3,$4,$1,$5,$6,$7,$8,$NF}' >> mm39_HoxD_custom_protein_coding.bed
# Make a version in which the HoxD cluster is represented instead of the individual Hoxd genes
hox_start=$(awk '$4 == "Hoxd13" {print $2}' mm39_HoxD_custom_protein_coding.bed)
hox_end=$(awk '$4 == "Hoxd1" {print $3}' mm39_HoxD_custom_protein_coding.bed)
awk -v OFS="\t" -v hox_start="$hox_start" -v hox_end="$hox_end" 'BEGIN {printed_hoxd = 0}{
    if (!($4 ~ /^Hoxd[0-9]+$/)){
        $9="255,255,255";print;
    } else if (printed_hoxd == 0) {
        print "chr2", hox_start, hox_end, "HoxD", "0", ".", hox_start, hox_end, "0,0,0";
        printed_hoxd = 1;
    }
}' mm39_HoxD_custom_protein_coding.bed > mm39_HoxD_custom_protein_coding_HoxDblack.bed
# Make a version of 3DOM regulated Hox genes
awk -v OFS="\t" '$4 ~ /^Hoxd([1-9]|10|11)$/ {$9="51,153,51";print}' mm39_HoxD_custom_protein_coding_details.bed > mm39_HoxD_custom_3DOMreg.bed
# Make a version of 5DOM regulated Hox genes
awk -v OFS="\t" '$4 ~ /^Hoxd1([0-3])$/ {$9="0,102,153";print}' mm39_HoxD_custom_protein_coding_details.bed > mm39_HoxD_custom_5DOMreg.bed

## Get element annotations from published primer sequences.
# If the primer sequences were not published, but the element was, we used the 20 bp flanking the element as an alternative and added a remark in the primer file.
# Format primer file to use in the perl script
awk -F',' 'NR>1 {OFS="\t"; print $1, $2, $3}' primers_transgenes_from_papers.txt > ${genomeDirectory}/primers_transgenes_from_papers.bed
# Get amplicon sequences
if [ ! -e "${genomeDirectory}/output_inSilico_PCR.txt" ]; then
    perl ../scripts/in_silico_PCR.pl -s ${genomeDirectory}/chr2.fa.gz -p ${genomeDirectory}/primers_transgenes_from_papers.bed > ${genomeDirectory}/output_inSilico_PCR.txt
else
    echo "${genomeDirectory}/output_inSilico_PCR.txt already exists, skipping Perl command."
fi
# Remove primer pairs which do not amplify
awk 'NR>1' ${genomeDirectory}/output_inSilico_PCR.txt | grep -v 'No amplification' > ${genomeDirectory}/output_inSilico_PCR_amp.txt
# Verify that each primer pair generates a unique amplicon
if awk 'NR>1' ${genomeDirectory}/output_inSilico_PCR_amp.txt | grep -v '_amp_1'; then
    echo "These primer pairs do not produce a unique amplicon."
else
    echo "Each primer pair produces a unique amplicon."
fi

# Get output in bed format.
awk -v OFS="\t" '{sub("_amp_1$", "", $1); print $2, $3, $3+$4, $1}' ${genomeDirectory}/output_inSilico_PCR_amp.txt > ${genomeDirectory}/mm39_HoxD_elements.bed
# Add lab annotations
cat annotations_HoxD_elements_mm39.bed >> ${genomeDirectory}/mm39_HoxD_elements.bed
# Add colors
awk -v OFS="\t" '{
  if ($3 < 74498654) {
    print $0, "0", ".",	$2, $3,	"0,102,153";
  } else {
    print $0, "0", ".",	$2, $3,	"51,153,51";
  }
}' ${genomeDirectory}/mm39_HoxD_elements.bed > mm39_HoxD_elements.bed

## Make different versions of lox annotations
awk -v OFS="\t" '($4 ~ /Itga6|nsi|attP/) {print $0}' annotations_lox_mm39.bed > mm39_HoxD_lox_Fig3.bed
awk -v OFS="\t" '($4 ~ /Rel|SB1|Atf2/) {print $0}' annotations_lox_mm39.bed > mm39_HoxD_lox_Fig4.bed

# Generate label tracks
for f in mm39_HoxD_custom_protein_coding_HoxDblack.bed mm39_HoxD_custom_protein_coding_details.bed; do 
    echo $f
    awk -v OFS="\t" '{print $1, $2, $2+1, $4}' $f > ${f/.bed/_start.bed}
done
