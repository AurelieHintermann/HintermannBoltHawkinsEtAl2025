library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggpubr)

# Get the raw data from the excel document
data <- read.xlsx("fluo_quantification/HCR_Fin_intensity_black_50.xlsx",
                  sep.names = " ")

# Normalize each value by the average in the wt
data <- data %>%
  mutate(relative_value = `Fin ROI Intensity` / mean(`Fin ROI Intensity`[Genotype == "wt"]))

# Use a better nomenclature for Genotype
pretty_genotype <- setNames(
  c("Wt", "Del(5DOM)"),
  c("wt", "del")
)

data$Genotype <- pretty_genotype[data$Genotype]

# Order the factors:
data$Genotype <- factor(data$Genotype, levels = unname(pretty_genotype))

# Compute the plot
ggplot(data, aes(x = Genotype, y = relative_value,
                 fill=Genotype)) +
  geom_boxplot(outliers = FALSE) + 
  ylim(0, 1.4) +
  geom_point(position=position_jitterdodge(0.2), alpha=1) +
  stat_compare_means(method = "t.test") +
  theme_classic(base_size = 16) +
  scale_fill_manual(values=c("white", "#E3D91B")) +
  ylab("Relative Intensity") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = NA, color = "black"))
ggsave("fluo_quantification/HCR_intensity.pdf",
       width = 3.5, height = 4.5)
