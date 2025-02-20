#! /bin/bash

if [ $0 != "-bash" ]; then
    source $(dirname $0)/../filePaths.sh
else
    echo "Source the file 'filePaths.sh' in the github directory"
fi

# This line is to adapt the conda to the shell
source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh

conda create -y -n HintermannBoltHawkinsEtAl2025 python=3.10.6 mamba
conda activate HintermannBoltHawkinsEtAl2025
mamba install --yes -c conda-forge -c bioconda --file ${gitHubDirectory}/HintermannBoltHawkinsEtAl2025.txt
conda deactivate
conda create -y -n HintermannBoltHawkinsEtAl2025_cufflinks python=3.6 mamba
conda activate HintermannBoltHawkinsEtAl2025_cufflinks
mamba install --yes -c conda-forge -c bioconda --file ${gitHubDirectory}/HintermannBoltHawkinsEtAl2025_cufflinks.txt
conda deactivate
