mkdir -p /scratch/ldelisle/maf
cd /scratch/ldelisle/maf
scriptDirectory=/scratch/ldelisle/HintermannBoltEtAl2024/maf_mm39_our_species/scripts
pathForTable=/scratch/ldelisle/HintermannBoltEtAl2024/maf_mm39_our_species/genomes_chr_select.txt

export PATH=$scriptDirectory:$PATH

while read l; do
    curgenome=$(echo $l | awk '{print $1}')
    chr=$(echo $l | awk '{print $2}')
    echo "Working on $curgenome"
    mkdir -p $curgenome
    if [ ! -e ${curgenome}/${curgenome}.fa.masked ]; then
        wget https://hgdownload.soe.ucsc.edu/goldenPath/${curgenome}/bigZips/${curgenome}.fa.masked.gz -P ${curgenome}/
        gunzip ${curgenome}/${curgenome}.fa.masked.gz
    fi
    if [ ! -e ${curgenome}/${curgenome}.2bit ]; then
        wget https://hgdownload.soe.ucsc.edu/goldenPath/${curgenome}/bigZips/${curgenome}.2bit -P ${curgenome}/
    fi
    faSize ${curgenome}/${curgenome}.fa.masked -detailed > ${curgenome}/${curgenome}.sizes
    echo "Splitting $curgenome"
    faSplit byName ${curgenome}/${curgenome}.fa.masked ${curgenome}/
    # To speed up we restrict to a single chromosome when we know
    if [ ! -z $chr ]; then
        mv ${curgenome}/${chr}.fa ${curgenome}/${chr}.fatokeep
        for f in $(find $curgenome -name "*.fa"); do
            rm $f
        done
        mv ${curgenome}/${chr}.fatokeep ${curgenome}/${chr}.fa
    fi
done < ${pathForTable}
