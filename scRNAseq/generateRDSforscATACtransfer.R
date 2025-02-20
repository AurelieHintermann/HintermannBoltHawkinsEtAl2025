options(scipen = 999)
options(stringsAsFactors = FALSE)
options(timeout = 10000)
rm(list = ls())
# Load packages for plot
if (!"ggplot2" %in% installed.packages()) {
  install.packages("ggplot2", repos = "https://stat.ethz.ch/CRAN/")
}
library(ggplot2)
if (!"ggpubr" %in% installed.packages()) {
  install.packages("ggpubr", repos = "https://stat.ethz.ch/CRAN/")
}
library(ggpubr)
if (!"dplyr" %in% installed.packages()) {
  install.packages("dplyr", repos = "https://stat.ethz.ch/CRAN/")
}
library(dplyr)
stages <- list(
  '14hpf' = 14
)
output.files <- paste0("scRNAseq/", names(stages), "_raw.rds")
names(output.files) <- names(stages)
if (any(!file.exists(output.files))) {
  # Load package required for single-cell
  if (!"Seurat" %in% installed.packages()) {
    install.packages("Seurat", repos = "https://stat.ethz.ch/CRAN/")
  }
  library(Seurat)
  # And package to read the supplementary table 2
  if (!"openxlsx" %in% installed.packages()) {
    install.packages("openxlsx", repos = "https://stat.ethz.ch/CRAN/")
  }
  library(openxlsx)
  # Download the data
  target.directory <- "scRNAseq/data/"
  dir.create(target.directory, showWarnings = FALSE, recursive = TRUE)
  url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE223922&format=file&file=GSE223922%5FSur2023%5Fcounts%2Emtx%2Egz"
  if (!file.exists(file.path(target.directory, "GSE223922_Sur2023_counts.mtx.gz"))) {
    download.file(url, file.path(target.directory, "GSE223922_Sur2023_counts.mtx.gz"))
  }
  url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE223922&format=file&file=GSE223922%5FSur2023%5Fcounts%5Fcols%5Fcells%2Etxt%2Egz"
  if (!file.exists(file.path(target.directory, "GSE223922_Sur2023_counts_cols_cells.txt.gz"))) {
    download.file(url, file.path(target.directory, "GSE223922_Sur2023_counts_cols_cells.txt.gz"))
  }
  url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE223922&format=file&file=GSE223922%5FSur2023%5Fcounts%5Frows%5Fgenes%2Etxt%2Egz"
  if (!file.exists(file.path(target.directory, "GSE223922_Sur2023_counts_rows_genes.txt.gz"))) {
    download.file(url, file.path(target.directory, "GSE223922_Sur2023_counts_rows_genes.txt.gz"))
  }
  url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE223922&format=file&file=GSE223922%5FSur2023%5Fmetadata%2Etsv%2Egz"
  if (!file.exists(file.path(target.directory, "GSE223922_Sur2023_metadata.tsv.gz"))) {
    download.file(url, file.path(target.directory, "GSE223922_Sur2023_metadata.tsv.gz"))
  }
  # Download the sup table 3
  url <- "https://www.cell.com/cms/10.1016/j.devcel.2023.11.001/attachment/98870784-3f9f-4028-99f9-7512b2bef58e/mmc3.xlsx"
  if (!file.exists(file.path(target.directory, "mmc3.xlsx"))) {
    stop("Please download mmc3.xlsx at ", url, " and save it to ", target.directory)
  }
  cluster.naming <- read.xlsx(xlsxFile = file.path(target.directory, "mmc3.xlsx"))
  cat("Building matrix..,")
  # Build matrix and make a Seurat object:
  expression_matrix <- ReadMtx(
    mtx = file.path(target.directory, "GSE223922_Sur2023_counts.mtx.gz"),
    features = file.path(target.directory, "GSE223922_Sur2023_counts_rows_genes.txt.gz"),
    cells = file.path(target.directory, "GSE223922_Sur2023_counts_cols_cells.txt.gz"),
    feature.column = 1
  )
  cat("Done\nBuilding Seurat Object")
  obj <- CreateSeuratObject(counts = expression_matrix)
  cat("Done\n")
  
  # Add all meta.data:
  meta.data.df <- read.delim(file.path(target.directory, "GSE223922_Sur2023_metadata.tsv.gz"))
  meta.data.df.enriched <- merge(meta.data.df, cluster.naming, by.x = "clust", by.y = "cluster", all.x = TRUE)
  rownames(meta.data.df.enriched) <- meta.data.df.enriched$cell
  obj <- AddMetaData(obj, meta.data.df.enriched)
  # Subset the obj to the selected stages
  for (my.label in names(stages)) {
    if (! file.exists(output.files[my.label])) {
      my.cells <- colnames(obj)[obj$stage.integer %in% stages[[my.label]]]
      cat("Keep", length(my.cells), "for", my.label, "\n")
      sub.obj <- subset(obj, cells = my.cells)
      saveRDS(sub.obj, output.files[my.label])
    }
  }
  with(
    subset(meta.data.df.enriched, stage.integer < 40), 
    table(stage.integer, clust == "endo.31")
  )
  # stage.integer FALSE  TRUE
  #          3    511     0
  #          4   2625     0
  #          5   5716     0
  #          6   1026     0
  #          7   4101     0
  #          8   6178     0
  #          9   5442     0
  #          10  7114     0
  #          11  1614     0
  #          12  4404     0
  #          14 12368    56
  #          16 11436    24
  #          18  8307     8
  #          21  9126     1
  #          24 19394     0
  #          26  2994     0
  #          28  5826     2
  #          30  5916     0
  #          32  7956     3
  #          34  8254     0
  #          36 30288     4
  #          38  4075     2
}


# Keep 12424 for 14hpf
