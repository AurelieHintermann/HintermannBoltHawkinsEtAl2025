options(scipen = 999)
options(stringsAsFactors = FALSE)
rm(list = ls())
# Load packages
if (!"reshape" %in% installed.packages()) {
  install.packages("reshape", repos = "https://stat.ethz.ch/CRAN/")
}
library(reshape)
if (!"ggplot2" %in% installed.packages()) {
  install.packages("ggplot2", repos = "https://stat.ethz.ch/CRAN/")
}
library(ggplot2)
if (!"pheatmap" %in% installed.packages()) {
  install.packages("pheatmap", repos = "https://stat.ethz.ch/CRAN/")
}
library(pheatmap)

# Set paths
file.with.FPKM <- commandArgs(TRUE)[1]
file.with.samples.plan <- commandArgs(TRUE)[2]
directory.with.plots <- commandArgs(TRUE)[3]
directory.with.data <- file.path(dirname(file.with.FPKM), "raw_values")
# Create non-existing output directories
if (!dir.exists(directory.with.plots)) {
  dir.create(directory.with.plots, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(directory.with.data)) {
  dir.create(directory.with.data, recursive = TRUE, showWarnings = FALSE)
}

FPKM.table <- read.delim(file.with.FPKM)
all.FPKM.values <- melt(FPKM.table)
all.FPKM.values$sample <- gsub("^FPKM_", "", all.FPKM.values$variable)

# Figure S8a:
# Get FPKM values for boxplot:
all.genes.needed <- paste0("Hoxd", 3:13)

FPKM.values <- subset(all.FPKM.values, gene_short_name %in% all.genes.needed)

# Get samplesplan
samplesplan <- read.delim(file.with.samples.plan)
samplesplan$genotype <- factor(samplesplan$genotype, levels = unique(samplesplan$genotype))

summary.df <- merge(FPKM.values, samplesplan)
summary.df$gene_short_name <- factor(summary.df$gene_short_name,
  levels = paste0("Hoxd", 13:3)
)

# Hoxd11 in Inv comes partially from a transgene
summary.df$alpha.value <- "show"
summary.df$alpha.value[summary.df$gene_short_name == "Hoxd11" & summary.df$genotype == "Inv"] <- "hide"

g <- ggplot(summary.df, aes(genotype, value, alpha = alpha.value)) +
  geom_boxplot(
    aes(fill = genotype, color = alpha.value),
    size = 0.2,
    outlier.shape = NA
  ) +
  geom_point(aes(shape = sex), size = 0.8) +
  facet_grid(. ~ gene_short_name, switch = "both") +
  theme_classic() +
  scale_fill_manual(
    values = c("wt" = "#b6bdc3", "Inv" = "#0e6dfc"),
    labels = c("wt" = "WT", "Inv" = "Inv(Itga6-nsi)d11lac")
  ) +
  scale_alpha_manual(values = c("hide" = 0.3, "show" = 1)) +
  scale_color_manual(values = c("hide" = "grey", "show" = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 21)) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(angle = 45, face = "italic", size = 7),
    element_line(linewidth = 0.1),
    axis.title = element_text(size = 7),
    # Adjust legend position
    legend.position = c(0.78, 0.85),
    legend.text = element_text(size = 7), # Adjust legend text size
    legend.key.size = unit(0.2, "cm"),
    legend.margin = margin(b = 0),
    strip.background = element_blank()
  ) +
  ylab("FPKM") +
  labs(shape = NULL, fill = NULL) +
  guides(colour = "none", alpha = "none")

ggsave(file.path(directory.with.plots, "FigS8a.pdf"), g, height = 7.1, width = 6.2, units = "cm")
# Export data:
write.table(cast(summary.df, sample + genotype + sex ~ gene_short_name, value = "value"),
  file.path(directory.with.data, paste0("FigS8a_raw_data.txt")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
