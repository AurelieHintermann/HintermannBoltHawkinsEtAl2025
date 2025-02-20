# Adapt code from https://github.com/Lan-lab/SPATAC-seq/blob/main/Zebrafish_Embryogenesis_single-cell_oPen_chromatin_Atlas_(ZEPA)/Figure%202.R
# Figure 2 
# cell type annotation 
# taking 14 hpf
dir.create("scATACseq/hpf14", showWarnings = FALSE)
setwd("scATACseq/hpf14/")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ArchR", quietly = TRUE)) devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
if (!requireNamespace("BSgenome.Drerio.UCSC.danRer11", quietly = TRUE)) BiocManager::install("BSgenome.Drerio.UCSC.danRer11")
if (!requireNamespace("org.Dr.eg.db", quietly = TRUE)) BiocManager::install("org.Dr.eg.db")
# To plot embedding:
if (!requireNamespace("hexbin", quietly = TRUE)) install.packages("hexbin")
#
###Zebrafish_
library(ArchR)
library(Seurat)
library("BSgenome.Drerio.UCSC.danRer11")
# library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(pheatmap)
library(ComplexHeatmap)
library("RColorBrewer")
library("circlize")
#
# loading annotation files
load("../danRer11_Lawson_genomeAnnotation_geneAnnotation.RData")
# Adapted:
chr_rm <- setdiff(
  seqlevels(geneAnnotation$genes),
  paste0("chr", 1:25)
)
##
# addArchRThreads(threads = 12) 
getArchRChrPrefix()
# load bed files
BedFiles <- c("../data/selected_14_sorted.bed.gz")
names(BedFiles) = c("14hpf")
# Create ArchR arrow file for downstream analysis
ArrowFiles <- createArrowFiles(
    inputFiles = BedFiles,
    sampleNames = names(BedFiles),
    minTSS = 4,
    minFrags = 500,
    offsetPlus = 0,
    offsetMinus = 0,
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation,
    addTileMat = T,
    addGeneScoreMat = T,
    excludeChr = chr_rm
)
# Found Gene Seqnames not in GenomeAnnotation chromSizes, Removing : chrM,chrUn_KN147632v2...
# Found Exon Seqnames not in GenomeAnnotation chromSizes, Removing : chrM,chrUn_KN147632v2,...
# Found TSS Seqnames not in GenomeAnnotation chromSizes, Removing : chrM,chrUn_KN147632v2,...
# ArchR logging to : ArchRLogs/ArchR-createArrows-31d9c047560313-Date-2024-08-09_Time-14-27-29.058385.log
# If there is an issue, please report to github with logFile!
# 2024-08-09 14:27:29.454639 : Batch Execution w/ safelapply!, 0 mins elapsed.
# Attempting to index ../data/selected_14_sorted.bed.gz as tabix..
# [W::tbx_parse1] Coordinate <= 0 detected. Did you forget to use the -0 option?
# [W::tbx_parse1] Coordinate <= 0 detected. Did you forget to use the -0 option?
# [W::tbx_parse1] Coordinate <= 0 detected. Did you forget to use the -0 option?
# [W::tbx_parse1] Coordinate <= 0 detected. Did you forget to use the -0 option?
# [W::tbx_parse1] Coordinate <= 0 detected. Did you forget to use the -0 option?
# [W::tbx_parse1] Coordinate <= 0 detected. Did you forget to use the -0 option?
# [W::tbx_parse1] Coordinate <= 0 detected. Did you forget to use the -0 option?
# [W::tbx_parse1] Coordinate <= 0 detected. Did you forget to use the -0 option?
# 2024-08-09 14:30:59.249352 : (14hpf : 1 of 1) Reading In Fragments from inputFiles (readMethod = tabix), 3.497 mins elapsed.
# 2024-08-09 14:31:00.664646 : (14hpf : 1 of 1) Tabix Bed To Temporary File, 3.52 mins elapsed.
# [W::tbx_parse1] Coordinate <= 0 detected. Did you forget to use the -0 option?
# [W::tbx_parse1] Coordinate <= 0 detected. Did you forget to use the -0 option?
# [W::tbx_parse1] Coordinate <= 0 detected. Did you forget to use the -0 option?
# [W::tbx_parse1] Coordinate <= 0 detected. Did you forget to use the -0 option?
# [W::tbx_parse1] Coordinate <= 0 detected. Did you forget to use the -0 option?
# [W::tbx_parse1] Coordinate <= 0 detected. Did you forget to use the -0 option?
# [W::tbx_parse1] Coordinate <= 0 detected. Did you forget to use the -0 option?
# [W::tbx_parse1] Coordinate <= 0 detected. Did you forget to use the -0 option?
# 2024-08-09 14:37:39.832743 : (14hpf : 1 of 1) Successful creation of Temporary File, 10.173 mins elapsed.
# 2024-08-09 14:37:39.836335 : (14hpf : 1 of 1) Creating ArrowFile From Temporary File, 10.173 mins elapsed.
# 2024-08-09 14:39:27.980318 : (14hpf : 1 of 1) Successful creation of Arrow File, 11.969 mins elapsed.
# 2024-08-09 14:42:06.157289 : (14hpf : 1 of 1) CellStats : Number of Cells Pass Filter = 47518 , 14.612 mins elapsed.
# 2024-08-09 14:42:06.159967 : (14hpf : 1 of 1) CellStats : Median Frags = 5005.5 , 14.612 mins elapsed.
# 2024-08-09 14:42:06.165502 : (14hpf : 1 of 1) CellStats : Median TSS Enrichment = 11.565 , 14.612 mins elapsed.
# 2024-08-09 14:42:22.569782 : (14hpf : 1 of 1) Adding Additional Feature Counts!, 14.885 mins elapsed.
# Error in order(typeScores, scores) : argument lengths differ
# In addition: Warning message:
# In completions$type == .rs.acCompletionTypes$DATAFRAME & !completions$context %in%  :
#   longer object length is not a multiple of shorter object length
# Error in order(typeScores, scores) : argument lengths differ
# In addition: Warning message:
# In completions$type == .rs.acCompletionTypes$DATAFRAME & !completions$context %in%  :
#   longer object length is not a multiple of shorter object length
# 2024-08-09 14:43:00.454805 : (14hpf : 1 of 1) Removing Fragments from Filtered Cells, 15.517 mins elapsed.
# 2024-08-09 14:43:04.389352 : (14hpf : 1 of 1) Creating Filtered Arrow File, 15.582 mins elapsed.
# 2024-08-09 14:47:41.216853 : (14hpf : 1 of 1) Finished Constructing Filtered Arrow File!, 20.196 mins elapsed.
# 2024-08-09 14:47:45.033275 : (14hpf : 1 of 1) Adding TileMatrix!, 20.26 mins elapsed.
# 2024-08-09 15:06:22.693305 : (14hpf : 1 of 1) Adding GeneScoreMatrix!, 38.887 mins elapsed.
# 2024-08-09 15:16:20.319966 : (14hpf : 1 of 1) Finished Creating Arrow File, 48.848 mins elapsed.
# ArchR logging successful to : ArchRLogs/ArchR-createArrows-31d9c047560313-Date-2024-08-09_Time-14-27-29.058385.log

# 
# I remove this because I think I don't need it:
# calculate Doublet score 
# doubScores <- addDoubletScores(
#   input = ArrowFiles,
#   k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
#   knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
#   LSIMethod = 1
# )
#

zhpf14 <- ArchRProject(ArrowFiles = ArrowFiles,
                       geneAnnotation = geneAnnotation,
                       genomeAnnotation = genomeAnnotation,
                       copyArrows = T)
#
zhpf14
# 47518 cells
#
df <- getCellColData(zhpf14, select = c("log10(nFrags)", "TSSEnrichment"))
df
p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 5, lty = "dashed") + geom_vline(xintercept = log10(1000), lty = "dashed")
p
# ggsave("uniqueFragsVSTSS.pdf", p, width = 6, height = 6)
#
idxPass <- which(zhpf14$TSSEnrichment >= 5 & zhpf14$nFrags >= 1000)
length(idxPass)
# All cells pass: 47518
cellsPass <- zhpf14$cellNames[idxPass]
zhpf14 <- zhpf14[cellsPass, ]
#
#
#
df <- getCellColData(zhpf14, select = c("log10(nFrags)", "TSSEnrichment"))
df
p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = log10(1000), lty = "dashed")
p
# Compute IterativeLSI with corCutOff of 0.5 (LM1) or 0.75 (LM1x)
zhpf14 <- addIterativeLSI(ArchRProj = zhpf14,
                          useMatrix = "TileMatrix", 
                          name = "IterativeLSI_Tile_LM1", 
                          iterations = 3, 
                          clusterParams = list( resolution = c(2,3), sampleCells = 10000, n.start = 10,maxClusters = 45), 
                          corCutOff = 0.5,
                          varFeatures = 50000, 
                          dimsToUse = 1:100,
                          LSIMethod = 1,
                          force=T)
# zhpf14 <- addIterativeLSI(ArchRProj = zhpf14,
#                           useMatrix = "TileMatrix", 
#                           name = "IterativeLSI_Tile_LM1x", 
#                           iterations = 3, 
#                           clusterParams = list( resolution = c(2,3), sampleCells = 10000, n.start = 10,maxClusters = 45), 
#                           corCutOff = 0.75,
#                           varFeatures = 50000, 
#                           dimsToUse = 1:100,
#                           LSIMethod = 1,
#                           force=T)
# save.image("zhpf14.RData")
barplot(zhpf14@reducedDims$IterativeLSI_Tile_LM1$corToDepth$scaled)
# barplot(zhpf14@reducedDims$IterativeLSI_Tile_LM1x$corToDepth$scaled)
#
zhpf14 <- addTSNE(ArchRProj = zhpf14, reducedDims = "IterativeLSI_Tile_LM1", name = "Tile_LM1_TSNE1", 
                  perplexity = 100, dimsToUse = 1:100,corCutOff = 0.4,maxIterations = 1000,force = T)
# zhpf14 <- addTSNE(ArchRProj = zhpf14, reducedDims = "IterativeLSI_Tile_LM1", name = "Tile_LM1_TSNE2", 
#                   perplexity = 100, dimsToUse = 1:100,corCutOff = 0.5,maxIterations = 1000,force = T)
# zhpf14 <- addTSNE(ArchRProj = zhpf14, reducedDims = "IterativeLSI_Tile_LM1", name = "Tile_LM1_TSNE3", 
#                   perplexity = 100, dimsToUse = 1:100,corCutOff = 0.8,maxIterations = 1000,force = T)
# save.image("zhpf14.RData")
#
zhpf14 <- addUMAP(ArchRProj = zhpf14, reducedDims = "IterativeLSI_Tile_LM1", name = "Tile_LM1_UMAP1", 
                  nNeighbors = 50, dimsToUse = 1:100, minDist = 0.3, metric = "cosine",force = T,corCutOff = 0.4)
# zhpf14 <- addUMAP(ArchRProj = zhpf14, reducedDims = "IterativeLSI_Tile_LM1", name = "Tile_LM1_UMAP2", 
#                   nNeighbors = 50, dimsToUse = 1:100, minDist = 0.3, metric = "cosine",force = T,corCutOff = 0.5)
# zhpf14 <- addUMAP(ArchRProj = zhpf14, reducedDims = "IterativeLSI_Tile_LM1", name = "Tile_LM1_UMAP3", 
#                   nNeighbors = 50, dimsToUse = 1:100, minDist = 0.3, metric = "cosine",force = T,corCutOff = 0.8)
#
# save.image("zhpf14.RData")
# zhpf14 <- addTSNE(ArchRProj = zhpf14, reducedDims = "IterativeLSI_Tile_LM1x", name = "Tile_LM1_TSNE1x", 
#                   perplexity = 100, dimsToUse = 1:100,corCutOff = 0.4,maxIterations = 1000,force = T)
# zhpf14 <- addTSNE(ArchRProj = zhpf14, reducedDims = "IterativeLSI_Tile_LM1x", name = "Tile_LM1_TSNE2x", 
#                   perplexity = 100, dimsToUse = 1:100,corCutOff = 0.5,maxIterations = 1000,force = T)
# zhpf14 <- addTSNE(ArchRProj = zhpf14, reducedDims = "IterativeLSI_Tile_LM1x", name = "Tile_LM1_TSNE3x", 
#                   perplexity = 100, dimsToUse = 1:100,corCutOff = 0.8,maxIterations = 1000,force = T)
# save.image("zhpf14.RData")
#
# zhpf14 <- addUMAP(ArchRProj = zhpf14, reducedDims = "IterativeLSI_Tile_LM1x", name = "Tile_LM1_UMAP1x", 
#                   nNeighbors = 50, dimsToUse = 1:100, minDist = 0.3, metric = "cosine",force = T,corCutOff = 0.4)
# zhpf14 <- addUMAP(ArchRProj = zhpf14, reducedDims = "IterativeLSI_Tile_LM1x", name = "Tile_LM1_UMAP2x", 
#                   nNeighbors = 50, dimsToUse = 1:100, minDist = 0.3, metric = "cosine",force = T,corCutOff = 0.5)
# zhpf14 <- addUMAP(ArchRProj = zhpf14, reducedDims = "IterativeLSI_Tile_LM1x", name = "Tile_LM1_UMAP3x", 
#                   nNeighbors = 50, dimsToUse = 1:100, minDist = 0.3, metric = "cosine",force = T,corCutOff = 0.8)
#
# save.image("zhpf14.RData")
#
plotEmbedding(ArchRProj = zhpf14, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_TSNE1")
# plotEmbedding(ArchRProj = zhpf14, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_TSNE2")
# plotEmbedding(ArchRProj = zhpf14, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_TSNE3")
plotEmbedding(ArchRProj = zhpf14, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_UMAP1")
# plotEmbedding(ArchRProj = zhpf14, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_UMAP2")
# plotEmbedding(ArchRProj = zhpf14, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_UMAP3")
# plotEmbedding(ArchRProj = zhpf14, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_TSNE1x")
# plotEmbedding(ArchRProj = zhpf14, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_TSNE2x")
# plotEmbedding(ArchRProj = zhpf14, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_TSNE3x")
# plotEmbedding(ArchRProj = zhpf14, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_UMAP1x")
# plotEmbedding(ArchRProj = zhpf14, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_UMAP2x")
# plotEmbedding(ArchRProj = zhpf14, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_UMAP3x")
#
#
#
zhpf14 <- addClusters(input = zhpf14,name = "Tile_LM1.C3",reducedDims = "IterativeLSI_Tile_LM1",corCutOff = 0.5,
                      resolution = 3,method = "Seurat", dimsToUse = 1:100,maxClusters= 300,force = T)
# zhpf14 <- addClusters(input = zhpf14,name = "Tile_LM1.C3x",reducedDims = "IterativeLSI_Tile_LM1x",corCutOff = 0.5,
#                       resolution = 3,method = "Seurat", dimsToUse = 1:100,maxClusters= 300,force = T)
# zhpf14 <- addClusters(input = zhpf14,name = "Tile_LM1.C5",reducedDims = "IterativeLSI_Tile_LM1",corCutOff = 0.5,
#                       resolution = 5,method = "Seurat", dimsToUse = 1:100,maxClusters= 300,force = T)
# zhpf14 <- addClusters(input = zhpf14,name = "Tile_LM1.C5x",reducedDims = "IterativeLSI_Tile_LM1x",corCutOff = 0.5,
#                       resolution = 5,method = "Seurat", dimsToUse = 1:100,maxClusters= 300,force = T)
#
#
# Added by me:
plotEmbedding(ArchRProj = zhpf14, colorBy = "cellColData", name = "Tile_LM1.C3", embedding = "Tile_LM1_TSNE1")
plotEmbedding(ArchRProj = zhpf14, colorBy = "cellColData", name = "Tile_LM1.C3", embedding = "Tile_LM1_UMAP1")
# plotEmbedding(ArchRProj = zhpf14, colorBy = "cellColData", name = "Tile_LM1.C3", embedding = "Tile_LM1_TSNE2")
# plotEmbedding(ArchRProj = zhpf14, colorBy = "cellColData", name = "Tile_LM1.C3x", embedding = "Tile_LM1_TSNE2x")
# plotEmbedding(ArchRProj = zhpf14, colorBy = "cellColData", name = "Tile_LM1.C5", embedding = "Tile_LM1_TSNE2")
# plotEmbedding(ArchRProj = zhpf14, colorBy = "cellColData", name = "Tile_LM1.C5x", embedding = "Tile_LM1_TSNE2x")
#
#
scRNAseq_14hpf <- readRDS("../../scRNAseq/14hpf_raw.rds")
scRNAseq_14hpf <- NormalizeData(scRNAseq_14hpf, normalization.method = "LogNormalize", scale.factor = 10000)
scRNAseq_14hpf <- FindVariableFeatures(scRNAseq_14hpf, selection.method = "vst", nfeatures = 3000)
scRNAseq_14hpf <- ScaleData(scRNAseq_14hpf, features = rownames(scRNAseq_14hpf))
scRNAseq_14hpf <- RunPCA(scRNAseq_14hpf, features = VariableFeatures(object = scRNAseq_14hpf))
scRNAseq_14hpf <- FindNeighbors(scRNAseq_14hpf, dims = 1:50)
scRNAseq_14hpf <- RunUMAP(scRNAseq_14hpf, dims = 1:50)
scRNAseq_14hpf <- RunTSNE(scRNAseq_14hpf, dims = 1:50)

scRNAseq_14hpf
#
# addArchRThreads(threads = 1) 
# This step took 90 minutes:
zhpf14 <- addGeneIntegrationMatrix(
  ArchRProj = zhpf14, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GIM_Tile_LSI.M1.1",
  reducedDims = "IterativeLSI_Tile_LM1",
  seRNA = scRNAseq_14hpf,
  addToArrow = T,
  force = TRUE,
  dimsToUse = 2:100,
  #corCutOff = 0.75,
  UMAPParams = list(n_neighbors = 50, min_dist = 0.3, metric = "cosine", verbose = FALSE),
  nGenes = 3000,
  sampleCellsATAC = length(zhpf14$cellNames),
  sampleCellsRNA = nrow(scRNAseq_14hpf[[]]),
  groupRNA = "clust",
  nameCell = "predictedCell_1",
  nameGroup = "predictedGroup_1",
  nameScore = "predictedScore_1"
)
#
#
# addArchRThreads(threads = 1) 
# hpf10.W
# zhpf14
# zhpf14 <- addGeneIntegrationMatrix(
#   ArchRProj = zhpf14, 
#   useMatrix = "GeneScoreMatrix",
#   matrixName = "GIM_Tile_LSI.M1.2",
#   reducedDims = "IterativeLSI_Tile_LM1",
#   seRNA = hpf10.W,
#   addToArrow = T,
#   force = TRUE,
#   dimsToUse = 2:100,
#   #corCutOff = 0.75,
#   UMAPParams = list(n_neighbors = 50, min_dist = 0.3, metric = "cosine", verbose = FALSE),
#   nGenes = 3000,
#   sampleCellsATAC = length(zhpf14$cellNames),
#   sampleCellsRNA = nrow(scRNAseq_14hpf[[]]),
#   groupRNA = "cell_type",
#   nameCell = "predictedCell_2",
#   nameGroup = "predictedGroup_2",
#   nameScore = "predictedScore_2")
#
#
#
#
plotEmbedding(ArchRProj = zhpf14, colorBy = "cellColData", name = "predictedGroup_1", embedding = "Tile_LM1_TSNE1")
plotEmbedding(ArchRProj = zhpf14, colorBy = "cellColData", name = "predictedGroup_1", embedding = "Tile_LM1_UMAP1")

# save.image("zhpf14.RData")
# Look at hoxd13a
plotEmbedding(ArchRProj = zhpf14, colorBy = "GeneScoreMatrix", name = "hoxd13a",
              embedding = "Tile_LM1_UMAP1",imputeWeights = getImputeWeights(zhpf14))


plotEmbedding(ArchRProj = zhpf14, colorBy = "GIM_Tile_LSI.M1.1", name = "hoxd13a",
              embedding = "Tile_LM1_UMAP1",imputeWeights = getImputeWeights(zhpf14))


df <- getCellColData(zhpf14, select = c("Tile_LM1.C3", "predictedGroup_1", "predictedScore_1"))
plot(density(df$predictedScore_1))
subset(df, predictedGroup_1 == "endo.31")
lines(density(df$predictedScore_1[df$predictedGroup_1 == "endo.31"]), col = "blue")

df2 <- getEmbedding(zhpf14, embedding = "Tile_LM1_UMAP1")
df <- cbind(df, df2[rownames(df), ])

df$tissue <- sapply(strsplit(df$predictedGroup_1, "\\."), "[[", 1)

ggplot(as.data.frame(df), aes(x = IterativeLSI_Tile_LM1.UMAP_Dimension_1, y = IterativeLSI_Tile_LM1.UMAP_Dimension_2)) +
  geom_point(aes(color = predictedGroup_1 == "endo.31")) +
  theme_minimal() +
  xlim(-10, 0) +
  ylim(-2, 2)

ggplot(as.data.frame(df), aes(x = IterativeLSI_Tile_LM1.UMAP_Dimension_1, y = IterativeLSI_Tile_LM1.UMAP_Dimension_2)) +
  geom_point(aes(color = predictedScore_1)) +
  theme_minimal() +
  xlim(-10, 0) +
  ylim(-2, 2)

ggplot(as.data.frame(df), aes(x = IterativeLSI_Tile_LM1.UMAP_Dimension_1, y = IterativeLSI_Tile_LM1.UMAP_Dimension_2)) +
  geom_point(aes(color = tissue)) +
  theme_minimal() +
  xlim(-10, 0) +
  ylim(-2, 2)


temp.df <- unique(scRNAseq_14hpf[[]][, c("clust", "identity.super")])
identity.super <- temp.df$identity.super
names(identity.super) <- temp.df$clust
zhpf14$identity.super <- identity.super[zhpf14$predictedGroup_1]

getGroupBW(
  ArchRProj = zhpf14,
  groupBy = "identity.super",
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)

saveRDS(zhpf14, "zhpf14.rds")

table(zhpf14$identity.super)

#         adaxial cells    cardiac mesenchyme                cloaca 
#                  1013                   904                    38 
#          diencephalon             epidermis           floor plate 
#                  1361                  4423                   529 
#        hatching gland        hemangioblasts             hindbrain 
#                   677                  1091                     1 
# intermediate mesoderm             ionocytes       mixed signature 
#                   477                   222                   556 
#         motor neurons          neural crest               neurons 
#                   901                  1296                   390 
#   non-neural ectoderm             notochord  olfactory epithelium 
#                    75                   859                    20 
#       optic primordia              periderm       pharyngeal arch 
#                  3100                   424                     7 
#      pharyngeal pouch     placodal ectoderm               placode 
#                   520                  2479                   449 
#   presomitic mesoderm           progenitors            sclerotome 
#                  2856                 12386                  1018 
#           slow muscle                somite      somitic mesoderm 
#                   233                   380                  2114 
#           spinal cord              tail bud         telencephalon 
#                  4378                   752                  1477 
#               unknown 
#                   112 
table(zhpf14$predictedGroup_1)

# axia.12 axia.13 axia.15 endo.27 endo.31  endo.7 epid.15 epid.22 epid.25 epid.26 
#     721     138     677     254      38     112    1476     488     520    1584 
#  epid.7  eye.26 glia.15  glia.2 glia.22  glia.6  glia.8 hema.19 hema.23 iono.21 
#    2839    3100     527       1       1     901    4378     800     291     222 
# mese.11  mese.9 musc.17 musc.21 musc.24 musc.25 musc.28  musc.3 musc.31  musc.4 
#     903       1     233      34    1013     556     752    2856     984    2114 
#  musc.8 neur.20 neur.21 neur.26 neur.27  neur.3 neur.30 neur.44  neur.6  neur.8 
#     380    1477    1361      75       7       1     389    1296   12132       1 
# otic.10 peri.16  pron.4  tast.4  tast.6 
#     449     424     477      20     515 

########################################################
zhpf14 <- readRDS("zhpf14.rds")

table(zhpf14$Tile_LM1.C3[zhpf14$predictedGroup_1 == "endo.31"])
table(zhpf14$predictedGroup_1[zhpf14$Tile_LM1.C3 == "C55"])
# axia.13 endo.27 endo.31 epid.25 glia.15  glia.8 hema.19 musc.25 musc.28  musc.3 neur.21  pron.4 
#      19     120      38       1       6       8      22       2      18       8       1      18 

# Try to visualize projection

ccaumap <- readRDS("ArchROutput/RNAIntegration/GIM_Tile_LSI.M1.1/Save-Block1-JointCCA.rds")
ggplot(as.data.frame(ccaumap)[sample(nrow(ccaumap)), ], aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Assay), size = 0.1) +
  theme_minimal() +
  annotate("rect", xmin = 0, xmax = 5, ymin = 5, ymax = 7.5, fill = NA, color = "black") +
  guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave("ccaumap_ATAC_RNA_highlight_Cloaca_region.pdf", height = 7, width = 7)


ggplot(as.data.frame(ccaumap), aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Group == "endo.31"), size = 0.1) +
  theme_minimal()

ccaumap.ordered <- rbind(subset(ccaumap, Group != "endo.31"), subset(ccaumap, Group == "endo.31"))

ggplot(as.data.frame(ccaumap.ordered), aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Group == "endo.31"), size = 0.1) +
  theme_minimal()

selected <- subset(
  as.data.frame(ccaumap.ordered),
  UMAP1 < 5 & UMAP1 > 0 & UMAP2 > 5 & UMAP2 < 7.5
)

g1 <- ggplot(selected, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Assay), size = 0.1) +
  theme_minimal() +
  guides(colour = guide_legend(override.aes = list(size = 3)))

g2 <- ggplot(selected, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Group == "endo.31"), size = 0.1) +
  theme_minimal() +
  scale_color_manual("is Cloaca", values = c('TRUE' = "red", 'FALSE' = "grey")) +
  guides(colour = guide_legend(override.aes = list(size = 3)))

g3 <- ggplot(subset(selected, Assay == "RNA"), aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Group == "endo.31"), size = 0.1) +
  theme_minimal() +
  ggtitle("RNA") +
  scale_color_manual("is Cloaca", values = c('TRUE' = "red", 'FALSE' = "grey")) +
  guides(colour = guide_legend(override.aes = list(size = 3)))

g4 <- ggplot(subset(selected, Assay == "ATAC"), aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Group == "endo.31"), size = 0.1) +
  theme_minimal() +
  ggtitle("ATAC") +
  scale_color_manual("is Cloaca", values = c('TRUE' = "red", 'FALSE' = "grey")) +
  guides(colour = guide_legend(override.aes = list(size = 3)))

g1.xfixed <- ggplot_gtable(ggplot_build(g1))
g1.xfixed$widths <- ggplot_gtable(ggplot_build(g2))$widths
library(ggpubr)
ggarrange(g1.xfixed, g2, g3, g4, ncol = 2, nrow = 2)
ggsave("ccaumap_ATAC_RNA_zoom_Cloaca_region.pdf", height = 7, width = 7)

sort(table(subset(selected, Assay == "RNA")$Group))
ggplot(subset(selected, Assay == "RNA"), aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Group)) +
  theme_minimal() +
  ggtitle("RNA")

# Get original cluster
# And package to read the metadata
if (!"openxlsx" %in% installed.packages()) {
install.packages("openxlsx", repos = "https://stat.ethz.ch/CRAN/")
}
library(openxlsx)  # Download the metadata

target.directory <- "../data/"
dir.create(target.directory, showWarnings = FALSE, recursive = TRUE)
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE243256&format=file&file=GSE243256%5FZEPA%5Fmetadata%2Exlsx"
if (!file.exists(file.path(target.directory, "GSE243256_ZEPA_metadata.xlsx"))) {
    download.file(url, file.path(target.directory, "GSE243256_ZEPA_metadata.xlsx"))
}
cell.metadata <- read.xlsx(xlsxFile = file.path(target.directory, "GSE243256_ZEPA_metadata.xlsx"))

cell.meta.14 <- subset(cell.metadata, stage == "14hpf")
rownames(cell.meta.14) <- paste0("14hpf#", cell.meta.14$cell_barcodes)
all(zhpf14$cellNames %in% rownames(cell.meta.14))
zhpf14$cell_type <- cell.meta.14[zhpf14$cellNames, "cell_type"]

saveRDS(zhpf14, "zhpf14.rds")

table(zhpf14$cell_type[zhpf14$predictedGroup_1 == "endo.31"])
# hpf14:Gut 
#        38 
