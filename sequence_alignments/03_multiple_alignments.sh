#!/bin/bash
conda activate BoltEtAl2023
cd /Users/Hintermann/Desktop/Bolt2023_test

annotation_file=annotations/HoxD_Elements_mm39_colored_noStrand9col.bed

wget -nc -O Genomes/homology/chr2_mm39.maf.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm39/multiz35way/maf/chr2.maf.gz

gunzip -k Genomes/homology/chr2_mm39.maf.gz

echo "
[x-axis]
fontsize = 5

[ATAC_bw]
file = /Users/Hintermann/switchdrive/Bolt_2023/toGEO/ATAC/ATAC_mm_E185_MUGS_rep1_macs_likeATAC_norm1.bw
title = MUGS_ATAC
show_labels = false
height = 1
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = darkblue

[spacer]
height = 0.1

[ChIP_narrowPeak]
file = /Users/Hintermann/switchdrive/Bolt_2023/toGEO/ATAC/ATAC_mm_E185_MUGS_rep1.narrowPeak
file_type = bed
color = darkblue
border_color = darkblue
labels = false
display = collapsed
line_width = 0.2
height = 0.2

[spacer]
height = 0.1

[ChIP_bw]
file = /Users/Hintermann/switchdrive/Bolt_2023/toGEO/ChIP/ChIP_H3K27ac_mm_E185_MUGS_rep1.bw
title = MUGS_H3K27ac
show_labels = false
height = 1
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = darkgreen

[spacer]
height = 0.1

[ChIP_narrowPeak]
file = /Users/Hintermann/switchdrive/Bolt_2023/toGEO/ChIP/ChIP_H3K27ac_mm_E185_MUGS_rep1.narrowPeak
file_type = bed
color = darkgreen
border_color = darkgreen
labels = false
display = collapsed
line_width = 0.2
height = 0.2

[spacer]
height = 0.1

[maf]
file = Genomes/homology/chr2_mm39.maf
file_type = maf
reference = mm39
file_index = Genomes/homology/chr2_mm39.maf.index
color_identical = black
color_mismatch = grey
color_gap = lightgrey
species_order = monDom5 xenTro9 galGal6 danRer11 petMar3
species_labels = Opossum Frog Chicken Zebrafish Lamprey
species_order_only = true
height = 1

[spacer]
height = 0.1

[features]
file = annotations/HoxD_Elements_mm39_colored_noStrand9col.bed
display = collapsed
color = bed_rgb
line_width = 1
border_color = bed_rgb
labels = true
fontsize = 6
height = .2
" > plots/Fig5_MAF.ini

# write inifile for 
for element in islandE CsB GT2; do 
    echo $element
    awk -F "\t" -v OFS="\t" -v el=$element '$0 ~ el"[[:space:]]0"' $annotation_file > tmp.bed
    start=$(cut -f 2 tmp.bed)
    start=$((start - 1000))
    stop=$(cut -f 3 tmp.bed)
    stop=$((stop + 1000))
    echo "#pygenometracks --tracks plots/Fig5_MAF.ini -o plots/Fig5_MAF_$(cut -f 4 tmp.bed).pdf --region $(cut -f 1 tmp.bed):$start-$stop --dpi 150 --width 6 --fontSize 6" > tmp.txt
    cat tmp.txt plots/Fig5_MAF.ini > plots/tmp_Fig5_MAF.ini
    mv plots/tmp_Fig5_MAF.ini plots/Fig5_MAF.ini
    pygenometracks --tracks plots/Fig5_MAF.ini -o plots/Fig5_MAF_$(cut -f 4 tmp.bed).pdf --region $(cut -f 1 tmp.bed):$start-$stop --dpi 150 --width 6 --fontSize 6
done

# There are lss species in multi alignment with mm39 than when Chase did it in mm10.
# Maybe it makes more sense to only show whether it is conserved or not in zf, using our conserved sequences annotated from our blast.

echo "[x-axis]
fontsize = 6

[ATAC_bw]
file = /Users/Hintermann/switchdrive/Bolt_2023/toGEO/ATAC/ATAC_mm_E185_MUGS_rep1_macs_likeATAC_norm1.bw
title = MUGS_ATAC
show_labels = false
height = 1
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = darkblue

[spacer]
height = 0.1

[ChIP_narrowPeak]
file = /Users/Hintermann/switchdrive/Bolt_2023/toGEO/ATAC/ATAC_mm_E185_MUGS_rep1.narrowPeak
file_type = bed
color = darkblue
border_color = darkblue
labels = false
display = collapsed
line_width = 0.2
height = 0.2

[spacer]
height = 0.1

[ChIP_bw]
file = /Users/Hintermann/switchdrive/Bolt_2023/toGEO/ChIP/ChIP_H3K27ac_mm_E185_MUGS_rep1.bw
title = MUGS_H3K27ac
show_labels = false
height = 1
file_type = bigwig
summary_method = max
nans_to_zeros = true
color = darkgreen

[spacer]
height = 0.1

[ChIP_narrowPeak]
file = /Users/Hintermann/switchdrive/Bolt_2023/toGEO/ChIP/ChIP_H3K27ac_mm_E185_MUGS_rep1.narrowPeak
file_type = bed
color = darkgreen
border_color = darkgreen
labels = false
display = collapsed
line_width = 0.2
height = 0.2

[spacer]
height = 0.1

[conserved sequences]
file = Genomes/homology/mm_conserved_elements_mm_zf.bed
title = conserved in zf
display = collapsed
color = black
line_width = 0
labels = false
height = 0.2

[spacer]
height = 0.1

[features]
file = annotations/HoxD_Elements_mm39_colored_noStrand9col.bed
display = collapsed
color = bed_rgb
line_width = 1
border_color = bed_rgb
labels = true
fontsize = 6
height = 0.2
fontsize = 6
" > plots/Fig5_trangenes.ini

for element in islandE CsB GT2; do 
    echo $element
    awk -F "\t" -v OFS="\t" -v el=$element '$0 ~ el"[[:space:]]0"' $annotation_file > tmp.bed
    start=$(cut -f 2 tmp.bed)
    start=$((start - 1000))
    stop=$(cut -f 3 tmp.bed)
    stop=$((stop + 1000))
    echo "#pygenometracks --tracks plots/Fig5_trangenes.ini -o plots/Fig5_trangenes_$(cut -f 4 tmp.bed).pdf --region $(cut -f 1 tmp.bed):$start-$stop --dpi 150 --width 6 --fontSize 6" > tmp.txt
    cat tmp.txt plots/Fig5_trangenes.ini > plots/Fig5_trangenes_tmp.ini
    mv plots/Fig5_trangenes_tmp.ini plots/Fig5_trangenes.ini
    pygenometracks --tracks plots/Fig5_trangenes.ini -o plots/Fig5_trangenes_$(cut -f 4 tmp.bed).pdf --region $(cut -f 1 tmp.bed):$start-$stop --dpi 150 --width 6 --fontSize 6
done

rm tmp*