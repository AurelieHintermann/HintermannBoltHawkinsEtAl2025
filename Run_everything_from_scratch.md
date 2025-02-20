# How to run everything from scratch

## Clone and make dirs

Clone the git repository

```bash
cd /scratch/ldelisle/
git clone git@github.com:AurelieHintermann/HintermannBoltHawkinsEtAl2025.git
```

Create dirs to be able to put fastqs

```bash
mkdir -p /scratch/ldelisle/HintermannBoltHawkinsEtAl/ATAC/fastq
mkdir -p /scratch/ldelisle/HintermannBoltHawkinsEtAl/cHiC/fastq
mkdir -p /scratch/ldelisle/HintermannBoltHawkinsEtAl/ChIP/fastq
mkdir -p /scratch/ldelisle/HintermannBoltHawkinsEtAl/CnR/fastq
mkdir -p /scratch/ldelisle/HintermannBoltHawkinsEtAl/RNAseq/fastq
```

Copy the fastqs into the right folders

Create a dir for HiC and TAD calling
```bash
mkdir -p /scratch/ldelisle/HintermannBoltHawkinsEtAl/
mkdir -p /scratch/ldelisle/HintermannBoltHawkinsEtAl/TAD_calling
```

## Run what is in prepare:

First modify what is in [filePaths.sh](./filePaths.sh)

Then run outside of jobs the small scripts in prepare:

```bash
cd /scratch/ldelisle/HintermannBoltHawkinsEtAl2025
bash prepare/00_create_conda_env.sh
bash prepare/01_get_fasta.sh
bash prepare/02a_custom_annotations_mm39.sh
bash prepare/02b_custom_annotations_danRer11.sh
```

Launch the slurm jobs to build indexes

```bash
source filePaths.sh
sbatch prepare/03_bowtie2_index.sh $gitHubDirectory $genomeDirectory
sbatch prepare/04_star_index.sh $gitHubDirectory $genomeDirectory
```

## Reproduce analyses:

### Sequence alignments:

While indexes are running, I run the sequence_alignments which does not require indexes:
```bash
bash sequence_alignments/01_homology_pipeline.sh
# To run the R script I need to activate the conda environment
source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
conda activate HintermannBoltHawkinsEtAl2025
Rscript sequence_alignments/02_homology_plots.R ./
conda deactivate
```

### MAF generation (1):

In order to evaluate the conservation of CsB, GT2 and islandE, we generated maf files for the selected species with the last version of the genome available.
To reduce the need in calculation power we restricted the alignment to mm39 chr2 as target and for query, see [this file](./maf_mm39_our_species/genomes_chr_select.txt).


#### Get all scripts

```bash
bash maf_mm39_our_species/00_get_all_scripts.sh
```

#### Get all fasta

```bash
bash maf_mm39_our_species/01_get_fasta.sh
```

#### Run lastz

```bash
sbatch maf_mm39_our_species/02_lastz.sh $PWD/maf_mm39_our_species/scripts $PWD/maf_mm39_our_species/genomes.txt /scratch/ldelisle/maf/mm39/chr2.fa
```

#### Run psl to maf

```bash
sbatch --dependency=aftercorr:13882911 maf_mm39_our_species/03_chaining_netting_maffing.sh $PWD/maf_mm39_our_species/scripts $PWD/maf_mm39_our_species/genomes.txt /scratch/ldelisle/maf/mm39/chr2.fa
```

#### Copy result to git
```bash
cp -r /scratch/ldelisle/maf/maf_per_species maf_mm39_our_species/
```


### MAF generation (2):

Another analysis was computed with cactus to evaluate the homology with gar and with the species themselves (duplication in the same genome)

```bash
bash sequence_alignments/10_cactus.sh
```

### ATAC

Run the sbatch script
```bash
sbatch ATAC/01_ATACseq_peak_coverage.sh $gitHubDirectory $genomeDirectory
```

Once done copy to GEO:
```bash
mkdir -p ${GEODirectory}/ATAC
cp /scratch/ldelisle/HintermannBoltHawkinsEtAl/ATAC/toGEO/* ${GEODirectory}/ATAC/
```

Copy versions to gitdir
```bash
cp /scratch/ldelisle/HintermannBoltHawkinsEtAl/ATAC/slurm*_01.out ${gitHubDirectory}/ATAC/
```

### ChIP

Run the sbatch scripts
```bash
sbatch ChIP/01_ChIP_SR.sh $gitHubDirectory $genomeDirectory
sbatch ChIP/01_ChIP_PE.sh $gitHubDirectory $genomeDirectory
```

Once done copy to GEO:
```bash
mkdir -p ${GEODirectory}/ChIP
cp /scratch/ldelisle/HintermannBoltHawkinsEtAl/ChIP/toGEO/* ${GEODirectory}/ChIP/
```

Copy versions to gitdir
```bash
cp /scratch/ldelisle/HintermannBoltHawkinsEtAl/ChIP/slurm*_01.out ${gitHubDirectory}/ChIP/
```

Run the CTCF orientation scripts from a temp folder:
```bash
mkdir -p /scratch/ldelisle/HintermannBoltHawkinsEtAl/ChIP/CTCF_orientation
cd /scratch/ldelisle/HintermannBoltHawkinsEtAl/ChIP/CTCF_orientation
bash ${gitHubDirectory}/ChIP/02_CTCF_orientation_mm.sh
```
I upload the result of insulatordb as requested.

Same with zebrafish:
```bash
bash ${gitHubDirectory}/ChIP/03_CTCF_orientation_zf.sh
```
I upload the result of insulatordb as requested.

### CUTandRUN

Run the sbatch script
```bash
sbatch CUTandRUN/CUTandRUN.sh $gitHubDirectory $genomeDirectory
```

Once done copy to GEO:
```bash
mkdir -p ${GEODirectory}/CnR
cp /scratch/ldelisle/HintermannBoltHawkinsEtAl/toGEO/* ${GEODirectory}/CnR/
```

Copy versions to gitdir
```bash
cp /scratch/ldelisle/HintermannBoltHawkinsEtAl/CnR/slurm*_01.out ${gitHubDirectory}/CUTandRUN/
```

### HiCandcHiC

### cHi-C
Run the sbatch script
```bash
sbatch HiCandcHiC/01_cHiC.sh $gitHubDirectory $genomeDirectory
```

Once done copy to GEO the 10kb cool and the valid pairs:
```bash
mkdir -p ${GEODirectory}/cHiC
cp /scratch/ldelisle/HintermannBoltHawkinsEtAl/cHiC/toGEO/*10kb* ${GEODirectory}/cHiC/
cp /scratch/ldelisle/HintermannBoltHawkinsEtAl/cHiC/toGEO/*S.valid* ${GEODirectory}/cHiC/
```

### Hi-C
Run the sbatch scripts
```bash
sbatch HiCandcHiC/02_HiC.sh $gitHubDirectory $genomeDirectory
sbatch --dependency afterok:12689931 HiCandcHiC/02b_HiC_merge_seq.sh $gitHubDirectory $genomeDirectory
sbatch --dependency afterok:12878316 HiCandcHiC/03_combine_cool.sh $gitHubDirectory
```

Once done copy to GEO the 10kb cool and the valid pairs:
```bash
mkdir -p ${GEODirectory}/HiC
cp /scratch/ldelisle/HintermannBoltHawkinsEtAl/HiC/toGEO/*10kb* ${GEODirectory}/HiC/
cp /scratch/ldelisle/HintermannBoltHawkinsEtAl/HiC/toGEO/*validPair* ${GEODirectory}/HiC/
cp /scratch/ldelisle/HintermannBoltHawkinsEtAl/HiC/allFinalFiles/cool/*merge*10kb.cool ${GEODirectory}/HiC/
```

### TAD calling

Once both are done, we can run the TAD calling
```bash
sbatch HiCandcHiC/04_TADcalling.sh $gitHubDirectory $GEODirectory
```

### RNA-seq

Run the sbatch script
```bash
sbatch RNAseq/01_RNAseq_SR.sh $gitHubDirectory $genomeDirectory
sbatch --dependency afterok:12689932 RNAseq/02_runRscripts.sh $gitHubDirectory $genomeDirectory
```
Copy versions to gitdir
```bash
cp /scratch/ldelisle/HintermannBoltHawkinsEtAl/RNAseq/slurm*_01.out ${gitHubDirectory}/RNAseq/
```

### scRNA-seq

The R script [plot.R](./scRNAseq/plot.R) has been run with a different R version and with different packages. I was run on a RStudio server (close to the docker image lldelisle/verse_with_more_packages:4.3.0_0 on dockerhub).
The output is a [table](./scRNAseq/outputs/endoderm_infos.txt) with all info in order to replot the UMAP and the gene expression without loading the Seurat object.
However the current version of the script uses ggpubr which is not in the conda environment.

```bash
# Outside of conda environment:
Rscript scRNAseq/plot.R 
```

In order to project cloaca cells from scRNAseq to scATACseq, a RDS file is generated containing only the 14hpf cells:

```bash
# Outside of conda environment:
Rscript scRNAseq/generateRDSforscATACtransfer.R
```


### RT-qPCR

Here again the R script [RTqPCR_plots.R](./RTqPCR/RTqPCR_plots.R) has been run with a different R with different packages. The current version of the script uses ggh4x which is not part of the conda environment.

```bash
# Outside of conda environment:
Rscript RTqPCR/RTqPCR_plots.R ./
```

### Cloaca measures

Here again the R script [cloaca_measures_plot.R](./cloaca_measures/cloaca_measures_plot.R) has been run with a different R with different packages. The current version of the script uses ggsignif and openxlsx which are not part of the conda environment.

```bash
# Outside of conda environment:
Rscript cloaca_measures/cloaca_measures_plot.R
```

### scATAC-seq

The idea is to use public scATAC and public scRNAseq where the Cloaca cells have been annotated to get a pseudo-bulk ATAC-seq profile of Cloaca.

First the barcodes of cells of 14hpf from the scATACseq are extracted:
```bash
# Outside of conda environment:
Rscript scATACseq/01_select_barcodes_14hpf.R
```

Then the fragments corresponding to these cells are extracted from the full fragment file and bgzip:
```bash
# Outside of conda environment:
bash scATACseq/02_extract_fragments_14hpf.sh
```

Then the analysis have been done with ArchR. The first step is to create an annotation file that uses the same gtf as in the scRNAseq:
```bash
# Outside of conda environment:
Rscript scATACseq/03_generate_annotation_files.R
```

Then the scRNAseq cells are transfered on the scATACseq cells of 14hpf.
```bash
# Outside of conda environment:
Rscript scATACseq/04_transfer_data_14hpf.R
```

### pyGenomeTracks plots

```bash
bash plots/01_plots.sh
```
