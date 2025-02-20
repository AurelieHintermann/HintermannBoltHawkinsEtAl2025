scriptDirectory=/scratch/ldelisle/HintermannBoltHawkinsEtAl2025/maf_mm39_our_species/scripts
mkdir -p $scriptDirectory
# wget https://genomewiki.ucsc.edu/images/9/93/RunLastzChain_sh.txt -nc -O ${scriptDirectory}/RunLastzChain.sh
# wget https://genomewiki.ucsc.edu/images/8/80/ConstructLiftFile_pl.txt -nc -O ${scriptDirectory}/ConstructLiftFile.pl

rsync -aP hgdownload.soe.ucsc.edu::genome/admin/exe/linux.x86_64/ ${scriptDirectory}/
