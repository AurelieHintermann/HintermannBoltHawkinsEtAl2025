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

sample <- "endoderm"
if (!file.exists(paste0("scRNAseq/outputs/", sample, "_infos.txt"))) {
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
  if (!file.exists(paste0("scRNAseq/", sample, ".rds"))) {
    if (!file.exists(paste0("scRNAseq/", sample, "_raw.rds"))) {
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

      # Build matrix and make a Seurat object:
      expression_matrix <- ReadMtx(
        mtx = file.path(target.directory, "GSE223922_Sur2023_counts.mtx.gz"),
        features = file.path(target.directory, "GSE223922_Sur2023_counts_rows_genes.txt.gz"),
        cells = file.path(target.directory, "GSE223922_Sur2023_counts_cols_cells.txt.gz"),
        feature.column = 1
      )
      obj <- CreateSeuratObject(counts = expression_matrix)

      # Add all meta.data:
      meta.data.df <- read.delim("scRNAseq/data/GSE223922_Sur2023_metadata.tsv.gz", row.names = 1)
      obj <- AddMetaData(obj, meta.data.df)

      ## FROM https://github.com/farrelllab/2023_Sur/blob/main/01_Pre-processing/merging_mama_dropseq.R
      # # Mitochondrial and ribosomal genes:
      # mito.genes <- grep("^mt-", rownames(obj@assays$RNA@counts), value=T)
      # ribo.genes <- intersect(rownames(obj@assays$RNA@counts), c("rpl18a", "rps16", "rplp2l", "rps13", "rps17", "rpl34",
      #                                                            "rpl13", "rplp0", "rpl36a", "rpl12", "rpl7a",
      #                                                            "rpl19", "rps2", "rps15a", "rpl3", "rpl27", "rpl23", "rps11",
      #                                                            "rps27a", "rpl5b", "rplp2", "rps26l", "rps10",
      #                                                            "rpl5a", "rps28", "rps8a", "rpl7", "rpl37", "rpl24", "rpl9",
      #                                                            "rps3a", "rps6", "rpl8", "rpl31", "rpl18", "rps27.2", "rps19",
      #                                                            "rps9", "rpl28", "rps7", "rpl7l1", "rps29", "rpl6",
      #                                                            "rps8b", "rpl10a", "rpl13a",
      #                                                            "rpl39", "rpl26", "rps24", "rps3", "rpl4",
      #                                                            "rpl35a", "rpl38", "rplp1", "rps27.1", "rpl15", "rps18", "rpl30",
      #                                                            "rpl11", "rpl14", "rps5", "rps21", "rpl10", "rps26",
      #                                                            "rps12", "rpl37.1", "rpl35", "rpl17", "rpl23a", "rps14", "rpl29",
      #                                                            "rps15", "rpl22", "rps23", "rps25", "rpl21",
      #                                                            "rpl22l1", "rpl36", "rpl32", "rps27l"))
      # obj[["percent.mt"]] <- PercentageFeatureSet(obj, features = mito.genes)
      # obj[["percent.ribo"]] <- PercentageFeatureSet(obj, features = ribo.genes)

      # message(paste0(Sys.time(), ": Saving results"))
      # saveRDS(obj, file="/data/CSD/zfext/results/merged_mama_DS_integrated_12012021_seurat.rds")

      # # Normalize and regress
      # message(paste0(Sys.time(), ": Normalizing Data"))
      # obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
      #
      # # message(paste0(Sys.time(), ": Saving results"))
      # # saveRDS(obj, file="/data/CSD/zfext/results/merged_mama_DS_integrated_12012021_seurat.rds")
      #
      # message(paste0(Sys.time(), ": Regressing out mito/ribo"))
      # obj <- ScaleData(obj, vars.to.regress = c("percent.mt", "percent.ribo"))
      #
      # # message(paste0(Sys.time(), ": Saving results"))
      # # saveRDS(obj, file="/data/CSD/zfext/results/merged_mama_DS_integrated_12012021_seurat.rds")
      #
      # # Variable Feature Selection
      # message(paste0(Sys.time(), ": Variable Feature Selection"))
      # obj <- FindVariableFeatures(obj, nfeatures = 2000)
      # # Save object and variable genes for downstream DR + clustering?
      # # message(paste0(Sys.time(), ": Saving results"))
      # # write(obj@assays$RNA@var.features, file=paste0("/data/CSD/zfext/results/merged_mama_DS_integrated_12012021_vargenes.txt"))
      # # saveRDS(obj, file="/data/CSD/zfext/results/merged_mama_DS_integrated_12012021_seurat.rds")
      #
      # # Dimensionality Reduction: PCA
      # message(paste0(Sys.time(), ": PCA"))
      # obj <- RunPCA(obj, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
      # obj <- JackStraw(obj, assay="RNA", dims=100)
      # obj <- ScoreJackStraw(obj, dims = 1:100)
      #
      # if (any(obj@reductions$pca@jackstraw@overall.p.values[,2] > 1e-3)) {
      #   dims.use <- min(which(obj@reductions$pca@jackstraw@overall.p.values[,2] > 1e-3)) - 1
      # } else {
      #   dims.use <- 30
      # }
      #
      # # Plot variable feature selection and PCA selection
      # # pdf(file= paste0("/data/CSD/zfext/results/merged_mama_DS_integrated_12012021_var.pdf"), width=8, height=8)
      # plot(VariableFeaturePlot(obj) + ggplot2::ggtitle(sample))
      # ElbowPlot(obj, ndims=100) + ggplot2::ggtitle(sample) + ggplot2::geom_vline(xintercept=dims.use+0.5, color='red')
      # # dev.off()
      #
      #
      # # Save object - temporarily in case of failure
      # # message(paste0(Sys.time(), ": Saving - temp"))
      # # saveRDS(obj, file="/data/CSD/zfext/results/merged_mama_DS_integrated_12012021_seurat.rds")
      #
      # # Dimensionality Reduction: UMAP
      # message(paste0(Sys.time(), ": UMAP"))
      # #nn.use <- floor(sqrt(ncol(obj@assays$RNA@counts))/10)*10
      # nn.use <- 50
      # obj <- RunUMAP(obj, dims=1:dims.use, n.neighbors=nn.use)

      # Temp skip
      # saveRDS(obj, "scRNAseq/all_cells_not_regressed.rds")

      #################################################
      ############### ENDODERM ########################
      #################################################

      ## ADAPTED FROM https://github.com/farrelllab/2023_Sur/blob/main/02_Clustering/subsets_crop_seurat.r

      # message(paste0(Sys.time(), ": Retrieve cell lists for ", sample))
      cells.list <- rownames(obj@meta.data)[which(obj@meta.data$tissue.name == sample)]

      # message(paste0(Sys.time(), ": Writing cells for", sample))
      # write(cells.list, paste0(cells.path, sample, "_cell_ids.txt"))
      cells.length <- length(cells.list)

      message(paste0(Sys.time(), ": ", cells.length, " cells are present in ", sample))

      # message(paste0(Sys.time(), ": Align cells with that in mama"))
      # cell.length.init <- length(cells.list)
      # cells.keep <- intersect(cells.list, colnames(mama@assays$RNA@counts))
      # cell.length.intersect <- length(cells.keep)
      # message(paste0(Sys.time(), ": ", cell.length.intersect, " of ", cell.length.init, " cells in the list were present in this object."))


      # Get relevant counts and metadata
      # message(paste0(Sys.time(), ": Getting relevant counts and metadata"))
      # subset.counts <- mama@assays$RNA@counts[, cells.keep]
      # subset.meta <- mama@meta.data[cells.keep, ]
      #
      # message(paste0(Sys.time(), ": Create Seurat object"))
      subset.counts <- obj@assays$RNA@counts[, cells.list]
      subset.meta <- meta.data.df[cells.list, ]
      obj.subset <- CreateSeuratObject(counts = subset.counts, meta.data = subset.meta, min.features = 0, min.cells = 2)

      if (!file.exists(paste0("scRNAseq/", sample, "_raw.rds"))) {
        saveRDS(obj.subset, paste0("scRNAseq/", sample, "_raw.rds"))
      }

      rm(list = c("obj"))
      shh <- gc()
      #
      # # Write out sample count matrices in MTX format needed for gene module analysis
      # message(paste0(Sys.time(), ": Savings counts as MTX"))
      # save.path <- paste0(counts.path, sample, ".mtx")
      # Matrix::writeMM(obj = obj@assays$RNA@counts, file = save.path)
      # write(rownames(obj@assays$RNA@counts), file = paste0(save.path, ".rn"))
      # write(colnames(obj@assays$RNA@counts), file = paste0(save.path, ".cn"))

      obj <- obj.subset
    } else {
      obj <- readRDS(paste0("scRNAseq/", sample, "_raw.rds"))
    }

    # Mitochondrial and ribosomal genes:
    mito.genes <- grep("^mt-", rownames(obj@assays$RNA@counts), value = T)
    ribo.genes <- intersect(rownames(obj@assays$RNA@counts), c(
      "rpl18a", "rps16", "rplp2l", "rps13", "rps17", "rpl34",
      "rpl13", "rplp0", "rpl36a", "rpl12", "rpl7a",
      "rpl19", "rps2", "rps15a", "rpl3", "rpl27", "rpl23", "rps11",
      "rps27a", "rpl5b", "rplp2", "rps26l", "rps10",
      "rpl5a", "rps28", "rps8a", "rpl7", "rpl37", "rpl24", "rpl9",
      "rps3a", "rps6", "rpl8", "rpl31", "rpl18", "rps27.2", "rps19",
      "rps9", "rpl28", "rps7", "rpl7l1", "rps29", "rpl6",
      "rps8b", "rpl10a", "rpl13a",
      "rpl39", "rpl26", "rps24", "rps3", "rpl4",
      "rpl35a", "rpl38", "rplp1", "rps27.1", "rpl15", "rps18", "rpl30",
      "rpl11", "rpl14", "rps5", "rps21", "rpl10", "rps26",
      "rps12", "rpl37.1", "rpl35", "rpl17", "rpl23a", "rps14", "rpl29",
      "rps15", "rpl22", "rps23", "rps25", "rpl21",
      "rpl22l1", "rpl36", "rpl32", "rps27l"
    ))
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, features = mito.genes)
    obj[["percent.ribo"]] <- PercentageFeatureSet(obj, features = ribo.genes)

    # Normalize and regress
    message(paste0(Sys.time(), ": Normalizing Data"))
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

    ## Regress out mito and ribo
    message(paste0(Sys.time(), ": Regressing out mito/ribo"))
    obj <- ScaleData(obj, vars.to.regress = c("percent.mt", "percent.ribo"))

    # Variable Feature Selection
    message(paste0(Sys.time(), ": Variable Feature Selection"))
    obj <- FindVariableFeatures(obj, nfeatures = 2000)

    # # Save object and variable genes for downstream DR + clustering?
    # message(paste0(Sys.time(), ": Saving results"))
    # write(obj@assays$RNA@var.features, file=paste0(var.path, sample, "_vargenes.txt"))
    # saveRDS(obj, file=paste0(seurat.path, sample, "_seurat.rds"))

    # Dimensionality Reduction: PCA
    message(paste0(Sys.time(), ": PCA"))
    obj <- RunPCA(obj, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
    obj <- JackStraw(obj, assay = "RNA", dims = 100)
    obj <- ScoreJackStraw(obj, dims = 1:100)

    if (any(obj@reductions$pca@jackstraw@overall.p.values[, 2] > 1e-3)) {
      dims.use <- min(which(obj@reductions$pca@jackstraw@overall.p.values[, 2] > 1e-3)) - 1
    } else {
      dims.use <- 30
    }

    message(paste0(Sys.time(), ": Dimensions used equal to :", dims.use))


    # # Plot variable feature selection and PCA selection
    # pdf(file=paste0(var.path, sample, "_var.pdf"), width=8, height=8)
    # plot(VariableFeaturePlot(obj) + ggplot2::ggtitle(sample))
    # ElbowPlot(obj, ndims=100) + ggplot2::ggtitle(sample) + ggplot2::geom_vline(xintercept=dims.use+0.5, color='red')
    # dev.off()


    # # Save object - temporarily in case of failure
    # message(paste0(Sys.time(), ": Saving - temp"))
    # saveRDS(obj, file=paste0(seurat.path, sample, "_seurat.rds"))

    # Dimensionality Reduction: UMAP
    message(paste0(Sys.time(), ": UMAP"))
    # nn.use <- floor(sqrt(ncol(obj@assays$RNA@counts))/10)*10
    nn.use <- 50
    obj <- RunUMAP(obj, dims = 1:dims.use, n.neighbors = nn.use)

    # Save object - temporarily in case of failure
    message(paste0(Sys.time(), ": Saving - temp"))

    DimPlot(obj, group.by = "clust")
    obj$cloaca <- "NO"
    obj$cloaca[obj$clust == "endo.31"] <- "YES"

    saveRDS(obj, paste0("scRNAseq/", sample, ".rds"))
    sessionInfo()
    # R version 4.3.0 (2023-04-21)
    # Platform: x86_64-pc-linux-gnu (64-bit)
    # Running under: Ubuntu 20.04.6 LTS

    # Matrix products: default
    # BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
    # LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

    # locale:
    #  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8
    #  [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
    #  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C
    # [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

    # time zone: Europe/Zurich
    # tzcode source: system (glibc)

    # attached base packages:
    # [1] stats     graphics  grDevices utils     datasets  methods   base

    # other attached packages:
    # [1] ggpubr_0.6.0       ggplot2_3.4.2      SeuratObject_4.1.3 Seurat_4.3.0

    # loaded via a namespace (and not attached):
    #   [1] RColorBrewer_1.1-3     rstudioapi_0.14        jsonlite_1.8.4
    #   [4] magrittr_2.0.3         spatstat.utils_3.0-2   ggbeeswarm_0.7.2
    #   [7] farver_2.1.1           ragg_1.2.5             fs_1.6.2
    #  [10] vctrs_0.6.2            ROCR_1.0-11            Cairo_1.6-0
    #  [13] memoise_2.0.1          spatstat.explore_3.1-0 rstatix_0.7.2
    #  [16] htmltools_0.5.5        usethis_2.1.6          broom_1.0.4
    #  [19] sctransform_0.3.5      parallelly_1.35.0      KernSmooth_2.23-21
    #  [22] htmlwidgets_1.6.2      ica_1.0-3              plyr_1.8.8
    #  [25] plotly_4.10.1          zoo_1.8-12             cachem_1.0.8
    #  [28] igraph_1.4.2           mime_0.12              lifecycle_1.0.3
    #  [31] pkgconfig_2.0.3        Matrix_1.5-4           R6_2.5.1
    #  [34] fastmap_1.1.1          fitdistrplus_1.1-11    future_1.32.0
    #  [37] shiny_1.7.4            digest_0.6.31          colorspace_2.1-0
    #  [40] patchwork_1.1.2        ps_1.7.5               tensor_1.5
    #  [43] irlba_2.3.5.1          pkgload_1.3.2          textshaping_0.3.6
    #  [46] labeling_0.4.2         progressr_0.13.0       fansi_1.0.4
    #  [49] spatstat.sparse_3.0-1  httr_1.4.5             polyclip_1.10-4
    #  [52] abind_1.4-5            compiler_4.3.0         remotes_2.4.2
    #  [55] withr_2.5.0            backports_1.4.1        carData_3.0-5
    #  [58] pkgbuild_1.4.0         ggsignif_0.6.4         MASS_7.3-60
    #  [61] sessioninfo_1.2.2      tools_4.3.0            vipor_0.4.5
    #  [64] lmtest_0.9-40          beeswarm_0.4.0         httpuv_1.6.9
    #  [67] future.apply_1.10.0    goftest_1.2-3          glue_1.6.2
    #  [70] callr_3.7.3            nlme_3.1-162           promises_1.2.0.1
    #  [73] grid_4.3.0             Rtsne_0.16             cluster_2.1.4
    #  [76] reshape2_1.4.4         generics_0.1.3         gtable_0.3.3
    #  [79] spatstat.data_3.0-1    tidyr_1.3.0            data.table_1.14.8
    #  [82] car_3.1-2              sp_1.6-0               utf8_1.2.3
    #  [85] spatstat.geom_3.1-0    RcppAnnoy_0.0.20       ggrepel_0.9.3
    #  [88] RANN_2.6.1             pillar_1.9.0           stringr_1.5.0
    #  [91] later_1.3.1            splines_4.3.0          dplyr_1.1.2
    #  [94] lattice_0.21-8         survival_3.5-5         deldir_1.0-6
    #  [97] tidyselect_1.2.0       miniUI_0.1.1.1         pbapply_1.7-0
    # [100] gridExtra_2.3          scattermore_1.0        devtools_2.4.5
    # [103] matrixStats_1.0.0      stringi_1.7.12         lazyeval_0.2.2
    # [106] codetools_0.2-19       tibble_3.2.1           cli_3.6.1
    # [109] uwot_0.1.14            systemfonts_1.0.4      xtable_1.8-4
    # [112] reticulate_1.28        munsell_0.5.0          processx_3.8.1
    # [115] Rcpp_1.0.10            globals_0.16.2         spatstat.random_3.1-4
    # [118] png_0.1-8              ggrastr_1.0.2          parallel_4.3.0
    # [121] ellipsis_0.3.2         prettyunits_1.1.1      profvis_0.3.8
    # [124] urlchecker_1.0.1       listenv_0.9.0          viridisLite_0.4.2
    # [127] scales_1.2.1           ggridges_0.5.4         leiden_0.4.3
    # [130] purrr_1.0.1            crayon_1.5.2           rlang_1.1.1
    # [133] cowplot_1.1.1
  } else {
    obj <- readRDS(paste0("scRNAseq/", sample, ".rds"))
  }

  hox13 <- intersect(
    paste0("hox", rep(letters[1:4], each = 2), 13, rep(letters[1:2], 4)),
    rownames(obj)
  )
  FeaturePlot(obj, features = hox13)
  Idents(obj) <- "clust"
  VlnPlot(obj, features = hox13)
  DotPlot(obj, features = hox13)

  temp.df <- FetchData(obj, vars = c("UMAP_1", "UMAP_2", "clust", "stage.integer", "stage.group", hox13))

  # Download the sup table 3
  target.directory <- "scRNAseq/data/"
  dir.create(target.directory, showWarnings = FALSE, recursive = TRUE)
  url <- "https://www.cell.com/cms/10.1016/j.devcel.2023.11.001/attachment/98870784-3f9f-4028-99f9-7512b2bef58e/mmc3.xlsx"
  if (!file.exists(file.path(target.directory, "mmc3.xlsx"))) {
    stop("Please download mmc3.xlsx at ", url, " and save it to ", target.directory)
  }

  cluster.naming <- read.xlsx(xlsxFile = file.path(target.directory, "mmc3.xlsx"))

  temp.df <- merge(temp.df, cluster.naming, by.x = "clust", by.y = "cluster", all.x = TRUE)

  write.table(temp.df, paste0("scRNAseq/outputs/", sample, "_infos.txt"), sep = "\t")
} else {
  temp.df <- read.delim(paste0("scRNAseq/outputs/", sample, "_infos.txt"))
}


# Get the name of genes from hox13 group
hox13 <- grep("^hox[a-d]13[ab]", colnames(temp.df), value = TRUE)

# Highlight cloaca and intestine
temp.df$cluster.highlight <- ""
temp.df$cluster.highlight[temp.df$identity.super %in% c("cloaca", "intestine")] <-
  temp.df$identity.super[temp.df$identity.super %in% c("cloaca", "intestine")]

# For intestine, get identity sub
temp.df$subcluster.highlight <- temp.df$cluster.highlight
temp.df$subcluster.highlight[temp.df$identity.super == "intestine"] <-
  paste0("Intestine\n", temp.df$identity.sub[temp.df$identity.super == "intestine"])

# Define colors
subcluster.highlight.colors <- c("grey", "red", colorRampPalette(c("#9cff9c", "#002900"))(length(unique(temp.df$subcluster.highlight)) - 2))
names(subcluster.highlight.colors) <- sort(unique(temp.df$subcluster.highlight))

# Set the limits to zoom:
xlims <- c(0, 10)
ylims <- c(-2.5, 7.5)

# Store ggplots to list:
gg.list <- list()

# First, all cells from endoderm
# With a rectangle around the zoomed area
gg.list[["Zoomout"]] <- ggplot(temp.df, aes(UMAP_1, UMAP_2)) +
  geom_point(size = .1, aes(color = .data[["identity.super"]])) +
  theme_classic() +
  geom_rect(
    aes(
      xmin = xlims[1],
      xmax = xlims[2],
      ymin = ylims[1],
      ymax = ylims[2]
    ),
    fill = NA,
    color = "black",
    linewidth = 1
  ) +
  theme(
    plot.title = element_text(size = 9, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.5, "lines")
  ) +
  labs(color = "") +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  ggtitle("Endoderm cells")

# Subset to cells in the rectangle
temp.df.subset <- subset(
  temp.df,
  UMAP_1 <= xlims[2] & UMAP_1 >= xlims[1] &
    UMAP_2 <= ylims[2] & UMAP_2 >= ylims[1]
)

# Readjust the levels to decrease the number of items in legend
for (colname in c("identity.super", "subcluster.highlight", "stage.group")) {
  temp.df.subset[, colname] <- factor(temp.df.subset[, colname])
}

# Make a summary df to hightlight Cloaca and Lysosome-rich
delta <- 1
summary.df <- temp.df.subset %>%
  group_by(subcluster.highlight) %>%
  summarise(
    min_x = min(UMAP_1),
    max_x = max(UMAP_1),
    mean_y = mean(UMAP_2)
  ) %>%
  mutate(
    x = ifelse(subcluster.highlight == "cloaca",
      max_x,
      min_x
    ),
    x_base = ifelse(subcluster.highlight == "cloaca",
      max_x + delta,
      min_x - delta
    ),
    keep = ifelse(subcluster.highlight == "cloaca" | grepl("lysosome", subcluster.highlight),
      TRUE,
      FALSE
    ),
    short_name = ifelse(grepl("EC", as.character(subcluster.highlight)),
      paste0("Intestine ", gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.|\n", "", as.character(subcluster.highlight), perl = TRUE)),
      as.character(subcluster.highlight)
    )
  )

gg.list[["identity.highlight"]] <-
  ggplot(temp.df.subset, aes(UMAP_1, UMAP_2, color = subcluster.highlight)) +
  geom_point(size = .1) +
  geom_segment(
    data = subset(summary.df, keep),
    aes(
      x = x_base, y = mean_y + delta,
      xend = x,
      yend = mean_y
    ),
    arrow = arrow(length = unit(0.5, "cm")),
    show.legend = FALSE
  ) +
  geom_text(
    data = subset(summary.df, keep),
    aes(
      label = short_name,
      x = x_base,
      y = mean_y + delta,
    ),
    vjust = -0.2,
    show.legend = FALSE
  ) +
  theme_classic() +
  scale_color_manual(
    "",
    values = subcluster.highlight.colors
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2), nrow = 5)) +
  theme(
    plot.title = element_text(size = 9, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.5, "lines")
  ) +
  ggtitle("Precised identity")

gg.list[["stage.group"]] <- ggplot(temp.df.subset, aes(UMAP_1, UMAP_2)) +
  geom_point(size = .1, aes(color = .data[["stage.group"]])) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 9, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.5, "lines")
  ) +
  labs(color = "") +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  ggtitle("Stage [hpf]")

gg.list.gene <- list()
for (gene in hox13) {
  gg.list.gene[[gene]] <- ggplot(temp.df.subset[order(temp.df.subset[, gene]), ], aes(UMAP_1, UMAP_2)) +
    geom_point(size = .1, aes(color = .data[[gene]])) +
    theme_classic() +
    scale_colour_gradientn(
      colours = c("grey", "orange", "red"),
      limits = c(0, 2),
      oob = scales::squish
    ) +
    ggtitle(gene) +
    theme(
      plot.title = element_text(size = 9, face = "bold.italic"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7, face = "bold.italic")
    ) +
    guides(colour = guide_colorbar(barheight = .5))
}
ggarrange(
  ggarrange(gg.list[["Zoomout"]],
    ggparagraph(text = "   "),
    ggarrange(
      gg.list[["identity.highlight"]],
      gg.list[["stage.group"]],
      ncol = 1, heights = c(1, .8)
    ),
    nrow = 1, widths = c(1.2, .1, 1)
  ),
  ggarrange(plotlist = gg.list.gene, common.legend = TRUE),
  ncol = 1,
  heights = c(1.3, 1)
)

ggsave("scRNAseq/outputs/FigS10.pdf", width = 19, height = 25, units = "cm")
