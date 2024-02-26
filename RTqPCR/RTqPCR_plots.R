options(scipen = 999)
options(stringsAsFactors = FALSE)
rm(list = ls())

# Load packages
if (!"dplyr" %in% installed.packages()) {
  install.packages("dplyr", repos = "https://stat.ethz.ch/CRAN/")
}
library(dplyr)
if (!"ggplot2" %in% installed.packages()) {
  install.packages("ggplot2", repos = "https://stat.ethz.ch/CRAN/")
}
library(ggplot2)
if (!"ggh4x" %in% installed.packages()) {
  install.packages("ggh4x", repos = "https://stat.ethz.ch/CRAN/")
}
library(ggh4x)

setwd(commandArgs(TRUE)[1])

# Set genotype order
genotype.order <- c(
  "wt", "Del(SB1-Atf2)", "Del(Rel5-SB1)", "Del(Rel1-Rel5)", "Inv/Inv", "Nsi-Atf2/Δ", "Rel5-SB1/Δ", "SB1-Atf2/Δ"
)

# Set paths
bed.file.with.gene.color <- "RTqPCR/HoxD_Elements_mm10_colored.bed"
download.file("https://raw.githubusercontent.com/lldelisle/scriptsForBoltEtAl2022/9644d5c4e3ea10b5dfae670fc08e25ddf2bbc78b/annotations/HoxD_Elements_mm10_colored.bed", bed.file.with.gene.color)
directory.with.plots <- "RTqPCR/plots"
directory.with.data <- "RTqPCR/raw_values"
# Create non-existing output directories
if (!dir.exists(directory.with.plots)) {
  dir.create(directory.with.plots, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(directory.with.data)) {
  dir.create(directory.with.data, recursive = TRUE, showWarnings = FALSE)
}

# Process colors
gene.color.info <- read.delim(bed.file.with.gene.color, header = FALSE)
gene.colors <- rgb(matrix(as.numeric(unlist(strsplit(gene.color.info$V9, ","))), byrow = TRUE, ncol = 3) / 256)
names(gene.colors) <- gene.color.info$V4

# Compute summary
summary.df.M <- NULL

for (file.with.RTqPCR in list.files(directory.with.data, full.names = TRUE)) {
  print(file.with.RTqPCR)

  # Read data file:
  raw.data.df <- read.delim(file.with.RTqPCR)

  raw.data.df$AnimalID <- raw.data.df$Sample

  # Check Nan data in wild type
  nan.wt.df <- raw.data.df %>%
    filter(is.na(Cq) & Genotype == "wt")

  # Remove wt animals with any NA
  filtered.data.df <- raw.data.df %>%
    filter(!AnimalID %in% nan.wt.df$AnimalID)

  # Average Cq per AnimalID/Plate and Target
  mean.Cq.df <- filtered.data.df %>%
    group_by(Target, AnimalID, Plate) %>%
    summarise(
      MeanCq = mean(Cq, na.rm = TRUE),
      n = sum(!is.na(Cq)), sd = sd(Cq, na.rm = TRUE)
    ) %>%
    filter(n >= 2 & sd < 0.4) %>% # These thresholds are arbitrary
    ungroup() %>%
    mutate(animalplate = paste(AnimalID, Plate, sep = "__"))

  # Clean AnimalId/Plate without HouseKeepting gene
  nan.HK <- mean.Cq.df %>%
    group_by(animalplate) %>%
    summarise(n = sum(Target == "Tbp")) %>%
    filter(n == 0)

  # Add delta Cq from Tbp
  dCq.df <- mean.Cq.df %>%
    filter(!animalplate %in% nan.HK$animalplate) %>%
    group_by(animalplate) %>%
    mutate(
      dCq = MeanCq - MeanCq[Target == "Tbp"]
    ) %>%
    ungroup() %>%
    select(Target, AnimalID, dCq)

  # Get delta delta Cq from average of Wt animals
  relExp.df <- dCq.df %>%
    left_join(filtered.data.df %>%
      select(AnimalID, Genotype) %>%
      unique()) %>%
    group_by(Target) %>%
    mutate(
      ddCq = dCq - mean(dCq[Genotype == "wt"]),
      relativeExpression = 2^-ddCq,
      exp = 2^-dCq,
      norm_relativeExpression = exp / mean(exp[Genotype == "wt"])
    ) %>%
    group_by(Target) %>%
    mutate(
      ddCqS = dCq - mean(dCq[Genotype == "wt"]),
      relativeExpressionS = 2^-ddCqS,
      norm_relativeExpressionS = exp / mean(exp[Genotype == "wt"])
    )

  selected.relExp.df <- relExp.df %>%
    filter(Genotype %in% genotype.order & Target %in% names(gene.colors) & !is.na(relativeExpression))
  selected.relExp.df$Target <- factor(selected.relExp.df$Target, levels = names(gene.colors))
  selected.relExp.df$Genotype <- factor(selected.relExp.df$Genotype, levels = genotype.order)

  # Plot all targets (FigS8)
  if (file.with.RTqPCR == "RTqPCR/raw_values/RTqPCR_E185_MUGS_CDOMdels_raw.txt") {
    set.seed(1)
    ggplot(selected.relExp.df, aes(x = Target, y = relativeExpressionS, fill = Target)) +
      geom_boxplot(outlier.colour = NA) +
      geom_jitter() +
      facet_nested(. ~ Genotype) +
      scale_fill_manual("Gene", values = gene.colors) +
      xlab("") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 7), axis.text.y = element_text(size = 7), axis.title.y = element_text(size = 7), legend.position = "none", element_line(linewidth = 0.1)) +
      ylab("relative expression") +
      scale_y_continuous(expand = c(0, 0), limits = c(0, ceiling(max(selected.relExp.df$relativeExpressionS))))

    ggsave(paste0(directory.with.plots, "/FigS8b_RTqPCR.pdf"), height = 2, width = 5.44)

    # Write file with quantifications displayed on FigS8b
    sink(paste0(directory.with.plots, "/FigS8b_RTqPCR.txt"))
    print(relExp.df %>%
      filter(Target == "Hoxd13") %>%
      group_by(Genotype) %>%
      summarise(
        Hoxd13_average_exp = mean(relativeExpression, na.rm = TRUE)
      ))
    sink()
  }

  # Plot Hoxd13 only
  if (file.with.RTqPCR != "RTqPCR/raw_values/RTqPCR_E185_MUGS_CDOMdels_raw.txt") {
    set.seed(1)
    ggplot(filter(selected.relExp.df, Target == "Hoxd13"), aes(x = Genotype, y = relativeExpressionS, fill = Target)) +
      geom_boxplot(outlier.colour = NA) +
      geom_jitter() +
      scale_fill_manual("Gene", values = gene.colors) +
      xlab("") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 7), axis.text.y = element_text(size = 7), axis.title.y = element_text(size = 7), legend.position = "none", element_line(linewidth = 0.1)) +
      ylab("relative expression") +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 2))

    ggsave(paste0(directory.with.plots, "/Fig3d_", unlist(strsplit(file.with.RTqPCR, "_"))[5], ".pdf"), height = 2.4, width = 1.1)
  }
}
