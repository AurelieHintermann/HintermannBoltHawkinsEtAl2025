Sinteract -m 50G -t 6:00:00
gitHubDirectory=/scratch/ldelisle/HintermannBoltHawkinsEtAl2025
wd=/scratch/ldelisle/HintermannBoltHawkinsEtAl/scATAC/
mkdir -p $wd
cd $wd
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE243256&format=file&file=GSE243256%5FZEPA%2EAll%2Esample%2Ebed%2Egz" -O GSE243256_ZEPA.All.sample.bed.gz
gunzip -k GSE243256_ZEPA.All.sample.bed.gz

# Select cells which are 14
awk '
NR == FNR{
    # store cell ids
    cells[$1] = 0
}
NR != FNR{
    if ($4 in cells) {
        print
    }
}
' ${gitHubDirectory}/scATACseq/cell_barcodes_14.txt GSE243256_ZEPA.All.sample.bed > selected_14.txt

# Takes too long:
# sort -k1,1 -k2,3n -k4,4 selected_14.txt > selected_14_sorted.bed
cat selected_14.txt | awk '{
print $0 > "tmp_"$1".bed"
}'
for file in tmp_*.bed; do
    sort -k2,3n -k4,4 $file >> selected_14_sorted.bed
done
module load gcc samtools
bgzip -k selected_14_sorted.bed

rm tmp_*.bed

# Copy it back to git:
mkdir -p $gitHubDirectory/scATACseq/data
cp selected_14_sorted.bed.gz $gitHubDirectory/scATACseq/data/
