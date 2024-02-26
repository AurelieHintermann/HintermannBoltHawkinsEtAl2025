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
if (!"ggsignif" %in% installed.packages()) {
  install.packages("ggsignif", repos = "https://stat.ethz.ch/CRAN/")
}
library(ggsignif)
# And package to read the supplementary table 2
if (!"openxlsx" %in% installed.packages()) {
  install.packages("openxlsx", repos = "https://stat.ethz.ch/CRAN/")
}
library(openxlsx)

data.directory <- "cloaca_measures"
excel_file_path <- file.path(data.directory, "cloaca_measures.xlsx")
measures <- read.xlsx(excel_file_path)
measures$Genotype <- factor(gsub(" ", "\n", measures$Genotype), levels = c("Wt", "hoxa13\nKO"))

set.seed(1)
ggplot(measures[measures$Measure == "Length", ], aes(x = Genotype, y = Value)) +
  geom_boxplot(fill = "lightgrey", outlier.colour = NA, size = 0.1) +
  geom_jitter(color = "blue", size = .5) +
  ggtitle("Length") +
  xlab("") +
  ylab("") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    plot.title = element_text(size = 7),
    legend.position = "none",
    element_line(linewidth = 0.1)
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  stat_compare_means(
    aes(label = paste0("p = ", after_stat(p.format))),
    method = "t.test",
    size = 2,
    label.y = 20,
    label.x = 1.5,
    hjust = 0.5
  ) +
  scale_y_continuous(breaks = seq(0, 100, 25), labels = c(0, rep("", 3), 100))
ggsave(paste0(data.directory, "/FigS11h_length.pdf"), height = 1.8, width = 1)

set.seed(1)
ggplot(measures[measures$Measure == "Width", ], aes(x = Genotype, y = Value)) +
  geom_boxplot(fill = "lightgrey", outlier.colour = NA, size = 0.1) +
  geom_jitter(color = "blue", size = .5) +
  ggtitle("Width") +
  xlab("") +
  ylab("") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        plot.title = element_text(size = 7),
        legend.position = "none",
        element_line(linewidth = 0.1)) +
  stat_compare_means(
    aes(label = paste0("p = ", after_stat(p.format))),
    method = "t.test",
    size = 2,
    label.y = 10,
    label.x = 1.5,
    hjust = 0.5
  ) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10), labels = c(0, rep("", 4), 50))
ggsave(paste0(data.directory, "/FigS11h_width.pdf"), height = 1.8, width = 1)
