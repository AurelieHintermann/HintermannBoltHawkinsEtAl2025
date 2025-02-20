options(timeout = 10000)
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ArchR", quietly = TRUE)) devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
if (!requireNamespace("BSgenome.Drerio.UCSC.danRer11", quietly = TRUE)) BiocManager::install("BSgenome.Drerio.UCSC.danRer11")
if (!requireNamespace("org.Dr.eg.db", quietly = TRUE)) BiocManager::install("org.Dr.eg.db")

library(ArchR)
library("BSgenome.Drerio.UCSC.danRer11")
library(org.Dr.eg.db)
library(GenomicFeatures)
genomeAnnotation <-
    createGenomeAnnotation(genome = "BSgenome.Drerio.UCSC.danRer11")
if (!file.exists("annotations/danRer11_Lawson_v4.3.2.gtf")) {
    download.file(
        url = "https://www.umassmed.edu/globalassets/lawson-lab/downloadfiles/v4.3.2.gtf",
        destfile = "annotations/danRer11_Lawson_v4.3.2.gtf"
    )
}
# This looses the gene_name:
# txdb <- makeTxDbFromGFF(
#     file =  "annotations/danRer11_Lawson_v4.3.2.gtf",
#     organism = "Danio rerio"
# )
gtf.gr <- import.gff("annotations/danRer11_Lawson_v4.3.2.gtf")
gtf.gr$gene_id <- gtf.gr$gene_name
txdb <- makeTxDbFromGRanges(gtf.gr)
geneAnnotation <- createGeneAnnotation(
    TxDb = txdb,
    OrgDb = org.Dr.eg.db,
    annoStyle = "SYMBOL"
)
save(genomeAnnotation, geneAnnotation, file = "scATACseq/danRer11_Lawson_genomeAnnotation_geneAnnotation.RData")
