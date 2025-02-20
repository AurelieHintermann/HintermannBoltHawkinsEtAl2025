#! /bin/bash
# Using pygenometracks version 3.8 HintermannBoltHawkinsEtAl2025

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

pathWithAnnotations="${gitHubDirectory}/annotations"
pathWithHiC="${GEODirectory}/HiC"
pathWithcHiC="${GEODirectory}/cHiC"
pathWithTADs="${gitHubDirectory}/HiCandcHiC/TADcalling"
pathWithCR="${GEODirectory}/CnR"
pathWithATAC="${GEODirectory}/ATAC"
pathWithChIP="${GEODirectory}/ChIP"
pathWithMAF="${gitHubDirectory}/maf_mm39_our_species/maf_per_species"
pathWithCactus3="${gitHubDirectory}/sequence_alignments/cactus/3way_HoxD_5DOM"
pathWithCactus13="${gitHubDirectory}/sequence_alignments/cactus/13way"
plottingDirectory="${gitHubDirectory}/plots/outputs"

mkdir -p $plottingDirectory
cd $plottingDirectory

# Figure 1a
ini_file=${plottingDirectory}/Fig1a.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --region chr9:1650000-2400000 --dpi 150 --decreasingXAxis --trackLabelFraction 0 --fontSize 7 --plotWidth 7.1"
echo -e "\n\n## Plotting ${ini_file/.ini/}"

echo "# ${pGT_cmd}

[spacer]
height = 0.1

[5DOM]
file = ${pathWithAnnotations}/danRer11_hoxDa_5DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .1
fontsize = 7

[3DOM]
file = ${pathWithAnnotations}/danRer11_hoxDa_3DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .1
fontsize = 7
overlay_previous = share-y

[spacer]
height = 0.1

[hline]
color = black
line_width = .3
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.15

[genes]
file = ${pathWithAnnotations}/danRer11_hoxDa_custom_protein_coding_HoxDblack.bed
file_type = bed
display = collapsed
color = bed_rgb
border_color = black
labels = false
height = 0.15
line_width = .5
overlay_previous = share-y

[spacer]
height = 0.2

[del3DOM]
file = ${pathWithAnnotations}/danRer11_hoxDa_del3DOM.bed
display = collapsed
color = red
border_color = none
height = .1
fontsize = 7

[spacer]
height = 0.1
"> $ini_file
$pGT_cmd

# Figure 1b
ini_file=${plottingDirectory}/Fig1b.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --region chr9:1650000-2400000 --dpi 150 --decreasingXAxis --trackLabelFraction 0 --fontSize 7 --plotWidth 6.6"
echo -e "\n\n## Plotting ${ini_file/.ini/}"

echo "# ${pGT_cmd}

[spacer]
height = 0.1

[5DOM]
file = ${pathWithAnnotations}/danRer11_hoxDa_5DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .1
fontsize = 7

[3DOM]
file = ${pathWithAnnotations}/danRer11_hoxDa_3DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .1
fontsize = 7
overlay_previous = share-y

[spacer]
height = 0.1

[hline]
color = black
line_width = .3
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.15

[genes]
file = ${pathWithAnnotations}/danRer11_hoxDa_custom_protein_coding_HoxDblack.bed
file_type = bed
display = collapsed
color = bed_rgb
border_color = black
labels = false
height = 0.15
line_width = .5
overlay_previous = share-y

[spacer]
height = 0.2

[del3DOM]
file = ${pathWithAnnotations}/danRer11_hoxDa_del5DOM.bed
display = collapsed
color = red
border_color = none
height = .1
fontsize = 7

[spacer]
height = 0.1
"> $ini_file
$pGT_cmd

# Figure 3c - part 1
ini_file=${plottingDirectory}/Fig3c1.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --region chr2:71600000-71700000 --dpi 150 --trackLabelFraction 0 --fontSize 7 --plotWidth 0.75"
echo -e "\n\n## Plotting ${ini_file/.ini/}"

echo "# ${pGT_cmd}
# Width was calculated to be proportional to Fig3c2. Default width is 40 cm. Default track label fraction is 0.05. ((0.95*40)*Fig3c1_region_size)/Fig3c2_region_size = 3.5
# then, to match .ai file, use 0.75 --plotWidth

[spacer]
height = 0.1

[hline]
color = black
line_width = .1
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.2

[genes]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding.bed
file_type = bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.2
line_width = .5
overlay_previous = share-y

[features labels]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding_details_start.bed
display = collapsed
color = none
line_width = 0
labels = true
fontsize = 7

[spacer]
height = 0.1

[lox sites]
file = ${pathWithAnnotations}/mm39_HoxD_lox_Fig3.bed
type = vlines
color = red
line_width = 1.5
alpha = 1
line_style = "dotted"
"> $ini_file
$pGT_cmd

# Figure 3c - part 2
ini_file=${plottingDirectory}/Fig3c2.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --region chr2:73630000-74720000 --dpi 150 --trackLabelFraction 0 --fontSize 7 --plotWidth 7"
echo -e "\n\n## Plotting ${ini_file/.ini/}"

echo "# ${pGT_cmd}

[spacer]
height = 0.1

[hline]
color = black
line_width = .1
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.2

[genes]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding.bed
file_type = bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.2
line_width = .5
overlay_previous = share-y

[lox sites]
file = ${pathWithAnnotations}/mm39_HoxD_lox_Fig3.bed
type = vlines
color = red
line_width = 1.5
alpha = 1
line_style = "dotted"

[spacer]
height = 0.1

[5DOM]
file = ${pathWithAnnotations}/mm39_HoxD_5DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .06
fontsize = 7
"> $ini_file
$pGT_cmd

# Figure 4a
ini_file=${plottingDirectory}/Fig4a.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --region chr2:73640000-74630000 --dpi 150 --trackLabelFraction 0.3 --fontSize 9 --plotWidth 16.7"
echo -e "\n\n## Plotting ${ini_file/.ini/}"

echo "# ${pGT_cmd}

[spacer]
height = 0.1

[5DOM]
file = ${pathWithAnnotations}/mm39_HoxD_5DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .08
fontsize = 7

[spacer]

[lox]
file = ${pathWithAnnotations}/mm39_HoxD_lox_Fig4.bed
file_type = bed
height = .3
color = bed_rgb
display = collapsed
line_width = 1
border_color = bed_rgb
labels = false
fontsize = 10

[ATAC_bw_mm_E185_MUGS_rep1]
file = ${pathWithATAC}/ATAC_mm_E185_MUGS_rep1.bw
title = ATAC_mm_E185_MUGS_rep1
height = 1.12
min_value = 0
max_value = 3.9
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = darkblue

[spacer]
height = 0.1

[ChIPac_H3K27ac_bw_mm_E185_MUGS_rep1]
file = ${pathWithChIP}/ChIP_H3K27ac_mm_E185_MUGS_rep1.bw
title = ChIP_H3K27ac_mm_E185_MUGS_rep1
height = 1.12
min_value = 0
max_value = 1.5
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = darkgreen

[spacer]
height = 0.1

[ChIPac_H3K27me3_bw_mm_E185_MUGS_rep1]
file = ${pathWithChIP}/ChIP_H3K27me3_mm_E185_MUGS_rep1.bw
title = ChIP_H3K27me3_mm_E185_MUGS_rep1
height = 1.12
min_value = 0
max_value = 6.8
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = darkred

[hline]
color = black
line_width = .5
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.2

[genes]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding.bed
file_type = bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.2
line_width = 1
overlay_previous = share-y

[elements]
file = ${gitHubDirectory}/annotations/mm39_HoxD_elements.bed
display = collapsed
color = bed_rgb
border_color = bed_rgb
height = 0.15
labels = false
overlay_previous = share-y

[features labels]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding_HoxDblack_start.bed
display = interleaved
color = none
line_width = 0
labels = true
fontsize = 7

[spacer]
height = 0.1
"> $ini_file
$pGT_cmd

# Figure S1a - part 1
ini_file=${plottingDirectory}/FigS1a1.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --region chr2:73600000-75550000 --dpi 150 --trackLabelFraction 0 --fontSize 7 --plotWidth 7"
echo -e "\n\n## Plotting ${ini_file/.ini/}"

echo "# ${pGT_cmd}
[hline]
color = black
line_width = .5
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.4

[genes]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding_HoxDblack.bed
file_type = bed
display = collapsed
color = bed_rgb
border_color = black
labels = false
height = 0.4
line_width = 1
overlay_previous = share-y

[spacer]
height = 0.05

[5DOM]
file = ${pathWithAnnotations}/mm39_HoxD_5DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .2
fontsize = 7

[3DOM]
file = ${pathWithAnnotations}/mm39_HoxD_3DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .2
fontsize = 7
overlay_previous = share-y

" > $ini_file
$pGT_cmd

# Figure S1a - part 2
ini_file=${plottingDirectory}/FigS1a2.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --region chr2:74497653-74560504 --dpi 150 --trackLabelFraction 0 --fontSize 7 --plotWidth 7"
echo -e "\n\n## Plotting ${ini_file/.ini/}"

echo "# ${pGT_cmd}

[hline]
color = black
line_width = .3
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.2

[genes]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding_details.bed
file_type = bed
display = collapsed
color = black
border_color = black
color_utr = black
labels = false
height = 0.2
arrowhead_included = true
overlay_previous = share-y

[genes]
file = ${pathWithAnnotations}/mm39_HoxD_custom_3DOMreg.bed
file_type = bed
display = collapsed
color = bed_rgb
border_color = bed_rgb
color_utr = bed_rgb
labels = false
height = 0.2
arrowhead_included = true
overlay_previous = share-y

[features labels]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding_details_start.bed
display = interleaved
color = none
line_width = 0
labels = true
fontsize = 7
"> $ini_file
$pGT_cmd

# Figure S1a - part 3
ini_file=${plottingDirectory}/FigS1a3.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --region chr2:73600000-75550000 --dpi 150  --fontSize 7 --plotWidth 7.1"
echo -e "\n\n## Plotting ${ini_file/.ini/}"

echo "# ${pGT_cmd}

[hline]
color = black
line_width = .3
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.2

[genes]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding_HoxDblack.bed
file_type = bed
display = collapsed
color = bed_rgb
border_color = black
labels = false
height = 0.2
line_width = .5
overlay_previous = share-y

[spacer]
height = 0.1

[5DOM]
file = ${pathWithAnnotations}/mm39_HoxD_5DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .1
fontsize = 7

[3DOM]
file = ${pathWithAnnotations}/mm39_HoxD_3DOM.bed
display = collapsed
color = #a3a3a3
border_color = none
height = .1
fontsize = 7
overlay_previous = share-y
" > $ini_file
$pGT_cmd

# Figure S1a - part 4
ini_file=${plottingDirectory}/FigS1a4.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --region chr2:74497653-74560504 --dpi 150 --trackLabelFraction 0 --fontSize 7 --plotWidth 7.1"
echo -e "\n\n## Plotting ${ini_file/.ini/}"

echo "# ${pGT_cmd}

[hline]
color = black
line_width = .3
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.2

[genes]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding_details.bed
file_type = bed
display = collapsed
color = black
border_color = black
color_utr = black
labels = false
height = 0.2
arrowhead_included = true
overlay_previous = share-y

[genes]
file = ${pathWithAnnotations}/mm39_HoxD_custom_5DOMreg.bed
file_type = bed
display = collapsed
color = bed_rgb
border_color = bed_rgb
color_utr = bed_rgb
labels = false
height = 0.2
arrowhead_included = true
overlay_previous = share-y

[features labels]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding_details_start.bed
display = interleaved
color = none
line_width = 0
labels = true
fontsize = 7
"> $ini_file
$pGT_cmd

# Figure S1b - part 1
ini_file=${plottingDirectory}/FigS1b1.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --region chr9:1650000-2400000 --dpi 150 --decreasingXAxis --trackLabelFraction 0 --fontSize 7 --plotWidth 7.1"
echo -e "\n\n## Plotting ${ini_file/.ini/}"

echo "# ${pGT_cmd}

[hline]
color = black
line_width = .3
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.2

[genes]
file = ${pathWithAnnotations}/danRer11_hoxDa_custom_protein_coding_HoxDblack.bed
file_type = bed
display = collapsed
color = bed_rgb
border_color = black
labels = false
height = 0.2
line_width = .5
overlay_previous = share-y
"> $ini_file
$pGT_cmd

# Figure S1b - part 2
ini_file=${plottingDirectory}/FigS1b2.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --region chr9:1932000-1995000 --dpi 150 --decreasingXAxis --trackLabelFraction 0 --fontSize 7 --plotWidth 7.1"
echo -e "\n\n## Plotting ${ini_file/.ini/}"

echo "# ${pGT_cmd}

[hline]
color = black
line_width = .3
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.2

[genes]
file = ${pathWithAnnotations}/danRer11_hoxDa_custom_preax.bed
file_type = bed
display = collapsed
color = bed_rgb
border_color = bed_rgb
color_utr = bed_rgb
labels = false
height = 0.2
arrowhead_included = true
overlay_previous = share-y

[features labels]
file = ${pathWithAnnotations}/danRer11_hoxDa_custom_preax_start.bed
display = interleaved
color = none
line_width = 0
labels = true
fontsize = 7
"> $ini_file
$pGT_cmd

# Figure S1b - part 3
ini_file=${plottingDirectory}/FigS1b3.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --region chr9:1932000-1995000 --dpi 150 --decreasingXAxis --trackLabelFraction 0 --fontSize 7 --plotWidth 7.1"
echo -e "\n\n## Plotting ${ini_file/.ini/}"

echo "# ${pGT_cmd}

[hline]
color = black
line_width = .3
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.2

[genes]
file = ${pathWithAnnotations}/danRer11_hoxDa_custom_postax.bed
file_type = bed
display = collapsed
color = bed_rgb
border_color = bed_rgb
color_utr = bed_rgb
labels = false
height = 0.2
arrowhead_included = true
overlay_previous = share-y

[features labels]
file = ${pathWithAnnotations}/danRer11_hoxDa_custom_preax_start.bed
display = interleaved
color = none
line_width = 0
labels = true
fontsize = 7
"> $ini_file
$pGT_cmd

# Figure S2 - part 1
ini_file=${plottingDirectory}/FigS2_1.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --region chr2:73600000-75550000 --dpi 150 --fontSize 7 --plotWidth 17.3"
echo -e "\n\n## Plotting ${ini_file/.ini/}"

echo "# ${pGT_cmd}

[CHiC_E185_MUGS_10kb_iced_mm39]
file = ${pathWithcHiC}/CHiC_mm_E185_MUGS.10kb.cool
title = CHiC E18.5 UGS
depth = 1000000
transform = no
file_type = hic_matrix
max_value = .04
height = 4.45

[spacer]
height = 0.05

[CHiC_E185_MUGS_10kb_iced_mm39_domains]
file = ${pathWithTADs}/CHiC_mm_E185_MUGS.10kb_240kb_domains.bed
border_color = none
color = black
line_width = 1.5
labels = false
height = .25

[ChIP_CTCF_E105_PT_rep1]
file = ${pathWithChIP}/CTCF_ChIP_mm_E105_PT_rep1.bw
title = CTCF E10.5 PT
height = .75
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = black

[spacer]
height = 0.05

[ChIP_CTCF_E105_PT_rep1_extended_noEMBL_colored_peakScorePlus500]
file = ${gitHubDirectory}/ChIP/CTCF_orientation/CTCF_ChIP_mm_E105_PT_rep1_extended_noEMBL_colored_peakScorePlus500.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.25

[spacer]
height = 0.05

[hline]
color = black
line_width = .3
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.2

[features]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding.bed
display = collapsed
color = bed_rgb
line_width = 1
border_color = bed_rgb
labels = false
height = 0.2
overlay_previous = share-y

[features labels]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding_HoxDblack_start.bed
display = interleaved
color = none
line_width = 0
fontsize = 7
labels = true
"> $ini_file
$pGT_cmd

# Figure S2 - part 2
ini_file=${plottingDirectory}/FigS2_2.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --region chr9:1650000-2400000 --dpi 300 --decreasingXAxis --fontSize 7 --plotWidth 17.3"
echo -e "\n\n## Plotting ${ini_file/.ini/}"

echo "# ${pGT_cmd}

[HiC_hpf2448merge_WE_10kb_iced_danRer11noAlt]
file = ${pathWithHiC}/HiC_hpf2448merge_WE_10kb.cool
title = HiC 24/48hpf full embryos
depth = 500000
transform = no
file_type = hic_matrix
#show_masked_bins = true
min_value = 0
max_value = .06
height = 5.75

[spacer]
height = 0.05

[HiC_hpf2448merge_WE_10kb_iced_danRer11noAlt_HoxDa_domains]
file = ${pathWithTADs}/HiC_hpf2448merge_WE_10kb_70kb_domains.bed
border_color = none
color = black
line_width = 1.5
labels = false
height = .25

[ChIP_CTCF_48hpf_WE_rep2]
file = ${pathWithChIP}/Franke_CTCF_ChIP_48hpf_rep2.bw
title = CTCF E10.5 PT
height = .75
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = black

[spacer]
height = 0.05

[Franke_ChIP_CTCF_allReps_noEMBL_colored]
file = ${gitHubDirectory}/ChIP/CTCF_orientation/Franke_CTCF_ChIP_allReps_noEMBL_colored.bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.25

[spacer]
height = 0.05

[hline]
color = black
line_width = .3
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.2

[features]
file = ${pathWithAnnotations}/danRer11_hoxDa_custom_protein_coding.bed
display = collapsed
color = bed_rgb
line_width = 1
border_color = bed_rgb
labels = false
height = 0.2
overlay_previous = share-y

[features labels]
file = ${pathWithAnnotations}/danRer11_hoxDa_custom_protein_coding_HoxDblack_start.bed
display = interleaved
color = none
line_width = 0
fontsize = 7
labels = true
"> $ini_file
$pGT_cmd

# Figure S3a - part 1
ini_file=${plottingDirectory}/FigS3a1.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --region chr2:73600000-75550000 --dpi 150 --fontSize 7 --plotWidth 13.4"
echo -e "\n\n## Plotting ${ini_file/.ini/}"

echo "# ${pGT_cmd}

[5DOM]
file = ${pathWithAnnotations}/mm39_HoxD_5DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .1
fontsize = 7

[3DOM]
file = ${pathWithAnnotations}/mm39_HoxD_3DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .1
fontsize = 7
overlay_previous = share-y

[spacer]
height = 0.1

[hline]
color = black
line_width = .3
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.15

[genes]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding.bed
file_type = bed
display = collapsed
color = bed_rgb
border_color = bed_rgb
labels = false
height = 0.15
line_width = .5
overlay_previous = share-y

[elements]
file = ${pathWithAnnotations}/mm39_HoxD_elements.bed
display = collapsed
color = bed_rgb
border_color = bed_rgb
height = 0.15
labels = false
overlay_previous = share-y

[spacer]
height = 0.1

[features labels]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding_HoxDblack_start.bed
display = interleaved
color = none
line_width = 0
labels = true
fontsize = 7

[spacer]
height = 0.1

[conserved sequences]
file = ${gitHubDirectory}/sequence_alignments/conserved_elements/mm_conserved_elements_mm_zf.bed
title = Conserved sequences
display = collapsed
color = black
line_width = 0.1
labels = false
height = .3
"> $ini_file
$pGT_cmd

# Figure S3a - part 2
ini_file=${plottingDirectory}/FigS3a2.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --region chr9:1650000-2400000 --dpi 150 --decreasingXAxis --fontSize 7 --plotWidth 13.4"
echo -e "\n\n## Plotting ${ini_file/.ini/}"

echo "# ${pGT_cmd}

[conserved sequences]
file = ${gitHubDirectory}/sequence_alignments/conserved_elements/zf_conserved_elements_mm_zf.bed
title = Conserved sequences
display = collapsed
color = black
line_width = 0.1
labels = false
height = .3

[spacer]
height = 0.1

[5DOM]
file = ${pathWithAnnotations}/danRer11_hoxDa_5DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .1
fontsize = 7

[3DOM]
file = ${pathWithAnnotations}/danRer11_hoxDa_3DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .1
fontsize = 7
overlay_previous = share-y

[spacer]
height = 0.1

[hline]
color = black
line_width = .3
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.15

[genes]
file = ${pathWithAnnotations}/danRer11_hoxDa_custom_protein_coding.bed
file_type = bed
display = collapsed
color = bed_rgb
border_color = bed_rgb
labels = false
height = 0.15
line_width = .5
overlay_previous = share-y

[elements]
file = ${pathWithAnnotations}/danRer11_hoxDa_elements.bed
display = collapsed
color = bed_rgb
border_color = bed_rgb
height = 0.15
labels = false
overlay_previous = share-y

[spacer]
height = 0.1

[features labels]
file = ${pathWithAnnotations}/danRer11_hoxDa_custom_protein_coding_HoxDblack_start.bed
display = interleaved
color = none
line_width = 0
labels = true
fontsize = 7
"> $ini_file
$pGT_cmd

# Figure S4c
ini_file=${plottingDirectory}/FigS4c.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --region chr9:1650000-2400000 --dpi 150 --decreasingXAxis --trackLabelFraction .25 --fontSize 7 --plotWidth 12.3"
echo -e "\n\n## Plotting ${ini_file/.ini/}"

echo "# ${pGT_cmd}

[CnR_danRer_16hpf_PT_H3K27ac_rep1]
file = ${pathWithCR}/CnR_H3K27ac_danRer_16hpf_PT_rep1.bw
title = CUT&RUN H3K27ac 16 hpf PT
height = .88
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = darkgreen

[spacer]
height = 0.1

[CnR_danRer_16hpf_PT_H3K27me3_rep1]
file = ${pathWithCR}/CnR_H3K27me3_danRer_16hpf_PT_rep1.bw
title = CUT&RUN H3K27me3 16 hpf PT
height = .88
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = darkred
max_value = 20

[spacer]
height = 0.1

[hline]
color = black
line_width = .3
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.2

[genes]
file = ${pathWithAnnotations}/danRer11_hoxDa_custom_protein_coding.bed
file_type = bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.2
line_width = .5
overlay_previous = share-y

[features labels]
file = ${pathWithAnnotations}/danRer11_hoxDa_custom_protein_coding_HoxDblack_start.bed
display = interleaved
color = none
line_width = 0
labels = true
fontsize = 7

[spacer]
height = 0.1

[5DOM]
file = ${pathWithAnnotations}/danRer11_hoxDa_5DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .08
fontsize = 7

[3DOM]
file = ${pathWithAnnotations}/danRer11_hoxDa_3DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .08
fontsize = 7
overlay_previous = share-y

[spacer]
height = 0.1

[hline]
color = black
line_width = .3
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.28

[genes]
file = ${pathWithAnnotations}/danRer11_hoxDa_custom_protein_coding_HoxDblack.bed
file_type = bed
display = collapsed
color = bed_rgb
border_color = black
labels = false
height = 0.28
line_width = .5
overlay_previous = share-y

"> $ini_file
$pGT_cmd

# Figure S7c
ini_file=${plottingDirectory}/FigS7c.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --region chr2:73630000-74720000 --dpi 150 --trackLabelFraction 0 --fontSize 7 --plotWidth 7"
echo -e "\n\n## Plotting ${ini_file/.ini/}"

echo "# ${pGT_cmd}

[spacer]
height = 0.1

[hline]
color = black
line_width = .1
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.2

[genes]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding.bed
file_type = bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.2
line_width = .5
overlay_previous = share-y

[lox sites]
file = ${pathWithAnnotations}/mm39_HoxD_lox_Fig3.bed
type = vlines
color = red
line_width = 1.5
alpha = 1
line_style = "dotted"

[spacer]
height = 0.1

[5DOM]
file = ${pathWithAnnotations}/mm39_HoxD_5DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .06
fontsize = 7

[spacer]
height = 0.1

[5DOM]
file = ${pathWithAnnotations}/mm39_HoxD_del.bed
display = collapsed
color = red
border_color = none
height = .1
fontsize = 7
"> $ini_file
$pGT_cmd

# Figure S7d
ini_file=${plottingDirectory}/FigS7d.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --region chr2:73630000-74720000 --dpi 150 --trackLabelFraction 0 --fontSize 7 --plotWidth 7"
echo -e "\n\n## Plotting ${ini_file/.ini/}"

echo "# ${pGT_cmd}

[spacer]
height = 0.1

[hline]
color = black
line_width = .1
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.2

[genes]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding.bed
file_type = bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.2
line_width = .5
overlay_previous = share-y

[lox sites]
file = ${pathWithAnnotations}/mm39_HoxD_lox_Fig3.bed
type = vlines
color = red
line_width = 1.5
alpha = 1
line_style = "dotted"

[spacer]
height = 0.1

[5DOM]
file = ${pathWithAnnotations}/mm39_HoxD_5DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .06
fontsize = 7

[spacer]
height = 0.1

[5DOM]
file = ${pathWithAnnotations}/mm39_HoxD_BAC.bed
display = collapsed
color = red
border_color = none
height = .1
fontsize = 7
"> $ini_file
$pGT_cmd

# Figure S8b
ini_file=${plottingDirectory}/FigS8b.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/.pdf} --region chr2:73640000-74630000 --dpi 150 --trackLabelFraction 0 --fontSize 7 --plotWidth 13"
echo -e "\n\n## Plotting ${ini_file/.ini/}"

echo "# ${pGT_cmd}

[spacer]
height = 0.1

[5DOM]
file = ${pathWithAnnotations}/mm39_HoxD_5DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .08
fontsize = 9

[spacer]
height = .7

[hline]
color = black
line_width = .5
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.2

[genes]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding.bed
file_type = bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.2
line_width = 1
overlay_previous = share-y

[elements]
file = ${pathWithAnnotations}/mm39_HoxD_elements.bed
display = collapsed
color = bed_rgb
border_color = bed_rgb
height = 0.2
labels = false
overlay_previous = share-y

[features labels]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding_HoxDblack_start.bed
display = interleaved
color = none
line_width = 0
labels = true
fontsize = 7

[spacer]
height = 0.1

[lox]
file = ${pathWithAnnotations}/mm39_HoxD_lox_Fig4.bed
file_type = bed
height = .3
color = bed_rgb
display = collapsed
line_width = 1
border_color = bed_rgb
labels = true
fontsize = 7

[spacer]
height = 0.1
" > $ini_file
$pGT_cmd


# Figure S9
ini_file=${plottingDirectory}/FigS9.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o  ${ini_file/.ini/_CsB.pdf} --region chr2:74280540-74283626 --dpi 150 --trackLabelFraction 0.3 --fontSize 7 --plotWidth 4.9"
pGT_cmd2="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/_islandE.pdf} --region chr2:74263903-74267008 --dpi 150 --trackLabelFraction 0.3 --fontSize 7 --plotWidth 4.9"
pGT_cmd3="pyGenomeTracks --tracks ${ini_file} -o ${ini_file/.ini/_GT2.pdf} --region chr2:74235672-74240601 --dpi 150 --trackLabelFraction 0.3 --fontSize 7 --plotWidth 4.9"

echo "# ${pGT_cmd}
# ${pGT_cmd2}
# ${pGT_cmd3}

[x-axis]
fontsize = 6

[ATAC_bw]
file = ${pathWithATAC}/ATAC_mm_E185_MUGS_rep1.bw
title = MUGS_ATAC
show_labels = false
height = .75
min_value = 0
max_value = 1.1
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = darkblue
number_of_bins = 200

[spacer]
height = 0.1

[ChIP_bw]
file = ${pathWithChIP}/ChIP_H3K27ac_mm_E185_MUGS_rep1.bw
title = MUGS_H3K27ac
show_labels = false
height = .75
min_value = .1
max_value = 0.8
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = darkgreen
max_value = 0.8
number_of_bins = 200

[spacer]
height = 0.1

[elements]
file = ${gitHubDirectory}/annotations/annotations_transgenes_mm39.bed
display = collapsed
color = #006699
border_color = none
height = 0.15
labels = true
fontsize = 7

[spacer]
height = 0.1
"> $ini_file

species=('ornAna2' 'galGal6' 'anoCar2' 'xenTro10' 'latCha1' 'fr3' 'danRer11' 'petMar3')
names=('Platypus' 'Chicken' 'Lizard' 'X_tropicalis' 'Coelacanth' 'Fugu' 'Zebraﬁsh' 'Lamprey')
names=('Pl' 'Ch' 'Li' 'XT' 'Coe' 'Fu' 'Ze' 'La')

for i in "${!species[@]}"; do
    echo "[MAF$i]
file = ${pathWithMAF}/mm39.chr2.${species[$i]}.maf
file_type = maf
reference = mm39
color_identical = black
color_mismatch = grey
color_gap = lightgrey
species_order = ${species[$i]}
species_labels = ${names[$i]}
height = 0.2
line_width = .1
title = ${species[$i]}
" >> $ini_file
done

echo -e "\n\n## Plotting ${ini_file/.ini/_CsB}"
$pGT_cmd
echo -e "\n\n## Plotting ${ini_file/.ini/_islandE}"
$pGT_cmd2
echo -e "\n\n## Plotting ${ini_file/.ini/_GT2}"
$pGT_cmd3

# Figure S10a
ini_file=${plottingDirectory}/FigS10a.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o  ${ini_file/.ini/.png} --region chr9:1900000-2400000 --dpi 300 --decreasingXAxis --trackLabelFraction .3 --fontSize 7 --plotWidth 16.7"

echo "# ${pGT_cmd}

[spacer]
height = 0.1

[maf]
file = ${pathWithCactus13}/vertebrates_13way_projected_Dr_chr9.maf
height = 5
reference = Danio_rerio
color_identical = black
color_mismatch = grey
color_gap = white
species_order = Danio_rerio Takifugu_rubripes Amia_calva Lepisosteus_oculatus Mus_musculus Ornithorhynchus_anatinus Gallus_gallus Anolis_carolinensis Xenopus_tropicalis Latimeria_chalumnae Scyliorhinus_canicula Leucoraja_erinacea Petromyzon_marinus
species_labels = zebrafish fugu bowfin gar mouse platypus chicken lizard frog coelacanth catshark skate lamprey
line_width = .1
title = whole genome alignments

[spacer]
height = 0.3

[hline]
color = black
line_width = .3
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.2

[genes]
file = ${pathWithAnnotations}/danRer11_hoxDa_custom_postax.bed
file_type = bed
display = collapsed
color = black
color = black
border_color = none
labels = false
height = 0.2
line_width = .5
overlay_previous = share-y

[features labels]
file = ${pathWithAnnotations}/danRer11_hoxDa_custom_protein_coding_HoxDblack_start.bed
display = interleaved
color = none
line_width = 0
labels = true
fontsize = 7

[elements]
file = ${pathWithAnnotations}/danRer11_hoxDa_elements.bed
display = collapsed
color = bed_rgb
border_color = bed_rgb
height = 0.15
labels = false

[elements labels]
file = ${pathWithAnnotations}/danRer11_hoxDa_elements.bed
display = collapsed
color = none
line_width = 0
labels = true
fontsize = 5

[spacer]
height = 0.3

[spacer]
height = 0.05

[5DOM]
file = ${pathWithAnnotations}/danRer11_hoxDa_5DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .08
fontsize = 7

[spacer]
height = 0.2

[x-axis]
fontsize = 6
" > ${ini_file}
$pGT_cmd

# Figure S10b
ini_file=${plottingDirectory}/FigS10b.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o  ${ini_file/.ini/.png} --region chr2:73640000-74630000 --dpi 300 --trackLabelFraction 0.3 --fontSize 7 --plotWidth 16.7"

echo "# ${pGT_cmd}

[x-axis]
fontsize = 6

[spacer]
height = 0.1

[5DOM]
file = ${pathWithAnnotations}/mm39_HoxD_5DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .08
fontsize = 7

[spacer]
height = 0.1

[hline]
color = black
line_width = .5
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.2

[genes]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding.bed
file_type = bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.2
line_width = 1
overlay_previous = share-y

[features labels]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding_HoxDblack_start.bed
display = interleaved
color = none
line_width = 0
labels = true
fontsize = 7
[spacer]
height = 0.1

[elements]
file = ${pathWithAnnotations}/mm39_HoxD_elements.bed
display = interleaved
color = bed_rgb
border_color = bed_rgb
height = 0.3
labels = true
fontsize = 5

[spacer]
height = 0.1

[maf]
file = ${pathWithCactus13}/vertebrates_13way_projected_Mm_chr2.maf
height = 5
reference = Mus_musculus
color_identical = black
color_mismatch = grey
color_gap = white
species_order = Mus_musculus Ornithorhynchus_anatinus Gallus_gallus Anolis_carolinensis Xenopus_tropicalis Latimeria_chalumnae Danio_rerio Takifugu_rubripes Amia_calva Lepisosteus_oculatus Scyliorhinus_canicula Leucoraja_erinacea Petromyzon_marinus
species_labels = mouse platypus chicken lizard frog coelacanth zebrafish fugu bowfin gar catshark skate lamprey
line_width = .1
title = whole genome alignments
"> $ini_file
$pGT_cmd

# Figure S10c
ini_file=${plottingDirectory}/FigS10c.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o  ${ini_file/.ini/.png} --region chr2:73640000-74630000 --dpi 300 --trackLabelFraction 0.3 --fontSize 7 --plotWidth 16.7"

echo "# ${pGT_cmd}

[x-axis]
fontsize = 6

[spacer]
height = 0.1

[5DOM]
file = ${pathWithAnnotations}/mm39_HoxD_5DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .08
fontsize = 7

[spacer]
height = 0.1

[hline]
color = black
line_width = .5
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.2

[genes]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding.bed
file_type = bed
display = collapsed
color = bed_rgb
border_color = none
labels = false
height = 0.2
line_width = 1
overlay_previous = share-y

[features labels]
file = ${pathWithAnnotations}/mm39_HoxD_custom_protein_coding_HoxDblack_start.bed
display = interleaved
color = none
line_width = 0
labels = true
fontsize = 7

[spacer]
height = 0.1

[elements]
file = ${pathWithAnnotations}/mm39_HoxD_elements.bed
display = interleaved
color = bed_rgb
border_color = bed_rgb
height = 0.3
labels = true
fontsize = 5

[spacer]
height = 0.1

[maf]
file = ${pathWithCactus3}/mm_zf_lo_projected_Mm.maf
height = 2
reference = Mus_musculus
line_width = 2.5
color_identical = black
color_mismatch = grey
color_gap = white
species_order = Mus_musculus Danio_rerio Lepisosteus_oculatus
species_labels = mouse zebrafish gar
title = HoxD 5DOM alignment
"> $ini_file
$pGT_cmd

# Compare pectoral fin and tailbud mesoderm
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE243256&format=file&file=GSE243256%5Fhpf24%5FTailbud%5Fmesoderm%2Ebw" -nc -O "GSE243256_hpf24_Tailbud_mesoderm.bw" 
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE243256&format=file&file=GSE243256%5Fhpf72%5FPectoral%5Ffin%5Ffield%2Ebw" -nc -O "GSE243256_hpf72_Pectoral_fin_field.bw"
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE243256&format=file&file=GSE243256%5Fhpf48%5FPectoral%5Ffin%5Ffield%2Ebw" -nc -O "GSE243256_hpf48_Pectoral_fin_field.bw"

# Figure S11
ini_file=${plottingDirectory}/FigS11.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o  ${ini_file/.ini/.pdf} --region chr9:1900000-2400000 --dpi 150 --decreasingXAxis --trackLabelFraction .25 --fontSize 7 --plotWidth 12.3"

echo "# ${pGT_cmd}

[x-axis]
fontsize = 7

[highlight names]
file = ${pathWithAnnotations}/danRer11_hoxDa_highlight_peaks.bed
color = none
height = 1.5
fontsize = 7

[spacer]
height = 0.1

[ATAC_danRer_30hpf_Head_rep1]
file = ${pathWithATAC}/ATAC_danRer_30hpf_Head_rep1.bw
title = ATAC Head 30hpf
height = .88
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = darkgrey
min_value = 0 

[spacer]
height = 0.1

[ATAC_danRer_30hpf_Cloaca_rep1]
file = ${pathWithATAC}/ATAC_danRer_30hpf_Cloaca_rep1.bw
title = ATAC Cloaca 30hpf
height = .88
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = brown
min_value = 0 

[spacer]
height = 0.1

[ATAC_danRer_30hpf_TailBud_rep1]
file = ${pathWithATAC}/ATAC_danRer_30hpf_TailBud_rep1.bw
title = ATAC TailBud 30hpf
height = .88
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = darkblue
min_value = 0 

[spacer]
height = 0.1

[hpf24_Tailbud_mesoderm]
file = GSE243256_hpf24_Tailbud_mesoderm.bw
title = scATAC Tailbud mesoderm 24hpf
height = .88
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = blue
min_value = 0 

[spacer]
height = 0.1

[scATAC_14hpf_endo.31]
file = ${gitHubDirectory}/scATACseq/hpf14/ArchROutput/GroupBigWigs/identity.super/cloaca-TileSize-100-normMethod-ReadsInTSS-ArchR.bw
title = scATAC Predicted cloaca 14hpf
height = .88
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = black
min_value = 0 

[spacer]
height = 0.1

[hpf48_Pectoral_fin_field]
file = GSE243256_hpf48_Pectoral_fin_field.bw
title = scATAC Pectoral fin 48hpf
height = .88
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = green

[spacer]
height = 0.1

[hpf72_Pectoral_fin_field]
file = GSE243256_hpf72_Pectoral_fin_field.bw
title = scATAC Pectoral fin 72hpf
height = .88
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = green
min_value = 0 

[spacer]
height = 0.1

[maf mouse]
file = ${gitHubDirectory}/maf_mm39_our_species/maf_per_species/mm39.chr2.danRer11.maf
file_index = ${gitHubDirectory}/maf_mm39_our_species/maf_per_species/mm39.chr2.danRer11.maf.indexDanRer
title = conservation with mouse
height = 0.2
reference = danRer11
color_identical = black
color_mismatch = grey
color_gap = lightgrey
line_width = .5
rasterize = true

[spacer]
height = 0.1

[hline]
color = black
line_width = .3
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.2
show_data_range = false

[genes]
file = ${pathWithAnnotations}/danRer11_hoxDa_custom_postax.bed
file_type = bed
display = collapsed
color = black
color = black
border_color = none
labels = false
height = 0.2
line_width = .5
overlay_previous = share-y

[features labels]
file = ${pathWithAnnotations}/danRer11_hoxDa_custom_protein_coding_HoxDblack_start.bed
display = interleaved
color = none
line_width = 0
labels = true
fontsize = 5

[elements]
file = ${pathWithAnnotations}/danRer11_hoxDa_elements.bed
display = collapsed
color = bed_rgb
border_color = bed_rgb
height = 0.15
labels = false

[elements labels]
file = ${pathWithAnnotations}/danRer11_hoxDa_elements.bed
display = collapsed
color = none
line_width = 0
labels = true
fontsize = 5

[spacer]
height = 0.05

[5DOM]
file = ${pathWithAnnotations}/danRer11_hoxDa_5DOM.bed
display = collapsed
color = bed_rgb
border_color = none
height = .08
fontsize = 7

[spacer]
height = 0.1

[del5DOM]
file = ${pathWithAnnotations}/danRer11_hoxDa_del5DOM.bed
display = collapsed
color = red
border_color = none
height = .1
fontsize = 7

[vhighlight]
file = ${pathWithAnnotations}/danRer11_hoxDa_highlight_peaks.bed
type = vhighlight
" > $ini_file
$pGT_cmd

# Figure S12a
ini_file=${plottingDirectory}/FigS12a.ini
pGT_cmd="pyGenomeTracks --tracks ${ini_file} -o  ${ini_file/.ini/.pdf} --region chr9:2068000-2080000 --dpi 150 --decreasingXAxis --trackLabelFraction .25 --fontSize 7 --plotWidth 12.3"

echo "# ${pGT_cmd}

[x-axis]
fontsize = 7

[spacer]
height = 0.1

[ATAC_danRer_30hpf_Head_rep1]
file = ${pathWithATAC}/ATAC_danRer_30hpf_Head_rep1.bw
title = ATAC Head 30hpf
height = .88
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = darkgrey
min_value = 0 

[spacer]
height = 0.1

[ATAC_danRer_30hpf_Cloaca_rep1]
file = ${pathWithATAC}/ATAC_danRer_30hpf_Cloaca_rep1.bw
title = ATAC Cloaca 30hpf
height = .88
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = brown
min_value = 0 

[spacer]
height = 0.1

[ATAC_danRer_30hpf_TailBud_rep1]
file = ${pathWithATAC}/ATAC_danRer_30hpf_TailBud_rep1.bw
title = ATAC TailBud 30hpf
height = .88
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = darkblue
min_value = 0 

[spacer]
height = 0.1

[hpf24_Tailbud_mesoderm]
file = GSE243256_hpf24_Tailbud_mesoderm.bw
title = scATAC Tailbud mesoderm 24hpf
height = .88
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = blue
min_value = 0 

[spacer]
height = 0.1

[scATAC_14hpf_endo.31]
file = ${gitHubDirectory}/scATACseq/hpf14/ArchROutput/GroupBigWigs/identity.super/cloaca-TileSize-100-normMethod-ReadsInTSS-ArchR.bw
title = scATACseq predicted cloaca 14hpf
height = .88
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = black
min_value = 0 

[spacer]
height = 0.1

[hpf48_Pectoral_fin_field]
file = GSE243256_hpf48_Pectoral_fin_field.bw
title = scATAC Pectoral fin 48hpf
height = .88
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = green

[spacer]
height = 0.1

[hpf72_Pectoral_fin_field]
file = GSE243256_hpf72_Pectoral_fin_field.bw
title = scATAC Pectoral fin 72hpf
height = .88
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = green
min_value = 0 

[spacer]
height = 0.1

[spacer]
height = 0.6

[maf]
file = ${gitHubDirectory}/maf_mm39_our_species/maf_per_species/mm39.chr2.danRer11.maf
height = 0.5
reference = danRer11
color_identical = black
color_mismatch = grey
color_gap = white
file_index = ${gitHubDirectory}/maf_mm39_our_species/maf_per_species/mm39.chr2.danRer11.maf.danRer.index
#species_order = Mus_musculus Danio_rerio
#species_labels = mouse zebrafish
line_width = .1
title = conservation

[spacer]
height = 0.3

[hline]
color = black
line_width = .3
min_value = 0
max_value = 1
file_type = hlines
y_values = .5
height = 0.2

[elements]
file = ${pathWithAnnotations}/danRer11_hoxDa_elements.bed
display = collapsed
color = bed_rgb
border_color = bed_rgb
height = 0.25
labels = true
title = liftOver
fontsize = 5
overlay_previous = share-y

[spacer]
height = 0.1

[dels]
file = ${pathWithAnnotations}/danRer11_delCsB.bed
title = CsB_guides
height = .3
color = red
border_color = red
fontsize = 7
display = collapsed
" > $ini_file
$pGT_cmd
