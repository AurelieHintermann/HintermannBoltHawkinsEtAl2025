# Bolt2023 Test

All scripts necessary to build figures from raw data in Hintermann, Bolt et al 2023. 

## To reproduce the analysis

Everything is described in detail in [this markdown](./Run_everyting_from_scratch.md).

Here is a summary:

### Setup

#### Conda environment

The analyses rely on a conda environment. To create it you need to:
- clone the git repository
- update the [filePaths.sh](./filePaths.sh)
- launch `bash prepare/00_create_conda_env.sh`.

#### Other requirements

Once this conda environment is created, you can run the bash scripts 01 and 02 in order to get fasta files. The 03 and 04 are designed to work with slurm scheduler and require two arguments: first, the gitDir (path where this repository has been cloned), then the genomeDirectory as you set it in the `filePaths.sh`.

### Preprocessing (fastq to processed files in GEO)

All slurm bash scripts used to generate the bigwig, narrowPeaks, counts, FPKM, cool files are available in the directory specific to each type of analysis: [ATAC](./ATAC/), [ChIP](./ChIP/), [CUTandRUN](./CUTandRUN/), [HiCandcHiC](./HiCandcHiC), [RNAseq](./RNAseq/).

### Plots

The script in [plots](./plots) uses pyGenomeTracks to generate final figures with coverages and heatmaps.

### Other analysis

Other analyses does not rely on fastqs: [sequence_alignments](./sequence_alignments/), [RTqPCR](./RTqPCR/), [scRNAseq](./scRNAseq/) and [cloaca_measures](./cloaca_measures/)

## Finally

All figures are symlinked into the [figures](./figures/) folder.
