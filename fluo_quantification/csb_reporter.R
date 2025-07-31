library(openxlsx)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggpubr)

# Get the raw data from the excel document
data_075 <- read.xlsx("fluo_quantification/CsB reporter intensity.xlsx",
                      "0.75 dpf")
data_075_long <- melt(data_075, id.vars = c("Fish", "Genotype"),
                      variable.name = "Tissue")
data_3 <- read.xlsx("fluo_quantification/CsB reporter intensity.xlsx",
                    "3 dpf", sep.names = " ")
data_3_long <- melt(data_3, id.vars = c("Fish", "Genotype"),
                    variable.name = "Tissue")
# Add the stage info
data_075_long$Stage <- "19 somites"
data_3_long$Stage <- "72 hpf"

# Merge both
data_all <- rbind(data_075_long, data_3_long)

# Normalize each Stage+Tissue value by the average in the Intact
data_all <- data_all %>%
  group_by(Stage, Tissue) %>%
  mutate(relative_value = value / mean(value[Genotype == "Intact"]))

# Order the factors:
data_all$Genotype <- factor(data_all$Genotype, levels = c("Intact", "Deleted"))
data_all$Tissue <- factor(data_all$Tissue,
                          levels = c("pectoral fin", "cloaca", "tail", "tailbud"))
# Compute the plot
ggplot(data_all, aes(x = Tissue, y = relative_value,
                     fill=Genotype)) +
  geom_boxplot() + 
  ylim(0.5, 2) +
  geom_point(position=position_jitterdodge(0.5), alpha=1) +
  facet_grid(. ~ Stage, scales = "free", space = "free", switch = "both") +
  stat_compare_means(aes(group = Genotype), label = "p.signif", method = "t.test",
                     label.y = 1.4, hide.ns = TRUE) +
  theme_classic(base_size = 16) +
  scale_fill_manual(values=c("white", "#F8766D"),
                    labels = paste("CsB", levels(data_all$Genotype))) +
  xlab("") +
  ylab("Relative Intensity") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.87, 0.8),
        legend.background = element_rect(fill = "white", color = "black"),
        # strip.background = element_blank(),
        strip.placement = "outside")
ggsave("fluo_quantification/CsB_reporter_intensity.pdf",
       width = 7, height = 4.5)

# Get the p-value values:
data_all$Stage_Tissue <- with(data_all, paste(Stage, Tissue, sep = "_"))
compare_means(value ~ Genotype, data_all, group.by = "Stage_Tissue", method = "t.test")
# Adding missing grouping variables: `Stage`, `Tissue`
# # A tibble: 5 × 9
#   Stage_Tissue        .y.   group1 group2       p p.adj p.format p.signif method
#   <chr>               <chr> <chr>  <chr>    <dbl> <dbl> <chr>    <chr>    <chr> 
# 1 19 somites_cloaca   value Intact Deleted 0.0221  0.11 0.022    *        T-test
# 2 19 somites_tailbud  value Intact Deleted 0.352   1    0.352    ns       T-test
# 3 72 hpf_pectoral fin value Intact Deleted 0.894   1    0.894    ns       T-test
# 4 72 hpf_cloaca       value Intact Deleted 0.274   1    0.274    ns       T-test
# 5 72 hpf_tail         value Intact Deleted 0.476   1    0.476    ns       T-test
