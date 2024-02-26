options(scipen = 999)
rm(list = ls())

library(ggplot2)
library(dplyr)

setwd(commandArgs(TRUE)[1])

conserved_sequence_file <- "sequence_alignments/conserved_elements/discontinuousMegablast_reciprocal_hits.txt"
conserved_sequence.pdf <- "sequence_alignments/outputs/FigS2b.pdf"
size_ratios.pdfs <- "sequence_alignments/outputs/FigS2c_"

dir.create(dirname(conserved_sequence.pdf), showWarnings = FALSE, recursive = TRUE)

# Get conserved sequences
cons.seq <- read.delim(conserved_sequence_file, header = TRUE)
# Get domain limits
# Transform all coordinate from bp to Mbp by /1e6
mm_5DOM_5 <- as.numeric(unlist(strsplit(readLines("annotations/mm39_HoxD_5DOM.bed"), "\t"))[2]) / 1e6
mm_5DOM_3 <- as.numeric(unlist(strsplit(readLines("annotations/mm39_HoxD_5DOM.bed"), "\t"))[3]) / 1e6
mm_3DOM_5 <- as.numeric(unlist(strsplit(readLines("annotations/mm39_HoxD_3DOM.bed"), "\t"))[2]) / 1e6
mm_3DOM_3 <- as.numeric(unlist(strsplit(readLines("annotations/mm39_HoxD_3DOM.bed"), "\t"))[3]) / 1e6
zf_5DOM_5 <- as.numeric(unlist(strsplit(readLines("annotations/danRer11_hoxDa_5DOM.bed"), "\t"))[3]) / 1e6
zf_5DOM_3 <- as.numeric(unlist(strsplit(readLines("annotations/danRer11_hoxDa_5DOM.bed"), "\t"))[2]) / 1e6
zf_3DOM_5 <- as.numeric(unlist(strsplit(readLines("annotations/danRer11_hoxDa_3DOM.bed"), "\t"))[3]) / 1e6
zf_3DOM_3 <- as.numeric(unlist(strsplit(readLines("annotations/danRer11_hoxDa_3DOM.bed"), "\t"))[2]) / 1e6

# For each conserved sequence, plot mouse coordinate on x and fish coordinate on y
synteny <- ggplot(cons.seq, aes(x = (as.numeric(mm_start) / 1e6), y = (as.numeric(-zf_start) / 1e6))) +
  geom_point(size = .8) +
  xlab("mouse coordinate [Mbp]") +
  ylab("zebrafish coordinate [Mbp]") +
  scale_y_continuous(labels = function(x) abs(x)) +
  theme_bw() +
  theme(element_line(linewidth = 0.1), axis.text = element_text(size = 7), axis.title = element_text(size = 7))

# Add lines at the position of domain limits
synteny <- synteny +
  geom_segment(aes(x = mm_5DOM_5, xend = mm_5DOM_3, y = -zf_5DOM_5, yend = -zf_5DOM_5), color = "#006699", linewidth = 0.1) +
  geom_segment(aes(x = mm_5DOM_3, xend = mm_5DOM_3, y = -zf_5DOM_5, yend = -zf_5DOM_3), color = "#006699", linewidth = 0.1) +
  geom_segment(aes(x = mm_5DOM_3, xend = mm_3DOM_5, y = -zf_5DOM_3, yend = -zf_5DOM_3), color = "#BE7FE6", linewidth = 0.1) +
  geom_segment(aes(x = mm_3DOM_5, xend = mm_3DOM_5, y = -zf_5DOM_3, yend = -zf_3DOM_5), color = "#BE7FE6", linewidth = 0.1) +
  geom_segment(aes(x = mm_3DOM_5, xend = mm_3DOM_3, y = -zf_3DOM_5, yend = -zf_3DOM_5), color = "#339933", linewidth = 0.1) +
  geom_segment(aes(x = mm_3DOM_3, xend = mm_3DOM_3, y = -zf_3DOM_5, yend = -zf_3DOM_3), color = "#339933", linewidth = 0.1)

ggsave(conserved_sequence.pdf, plot = synteny, units = "cm", width = 6.5, height = 4.8)

# Plot domain size ratios [Mbp]
assembly_mm39 <- 2.7 * 1e3
assembly_danRer11 <- 1.4 * 1e3

# Initiate data frame with assembly sizes
domain_sizes_Mbp <- data.frame(c("mm39", "assembly", assembly_mm39), c("danRer11", "assembly", assembly_danRer11), colnames = c("assembly", "domain", "size"))
domain_sizes_Mbp <- data.frame(c("mm39", "danRer11"), c("assembly", "assembly"), c(assembly_mm39, assembly_danRer11))
colnames(domain_sizes_Mbp) <- c("assembly", "domain", "size")

# Fill in data frame with sizes of the domains
domain_sizes_Mbp <- rbind(domain_sizes_Mbp, c("mm39", "5DOM", abs(mm_5DOM_5 - mm_5DOM_3)))
domain_sizes_Mbp <- rbind(domain_sizes_Mbp, c("danRer11", "5DOM", abs(zf_5DOM_5 - zf_5DOM_3)))
domain_sizes_Mbp <- rbind(domain_sizes_Mbp, c("mm39", "3DOM", abs(mm_3DOM_5 - mm_3DOM_3)))
domain_sizes_Mbp <- rbind(domain_sizes_Mbp, c("danRer11", "3DOM", abs(zf_3DOM_5 - zf_3DOM_3)))
domain_sizes_Mbp <- rbind(domain_sizes_Mbp, c("mm39", "cluster", abs(mm_5DOM_3 - mm_3DOM_5)))
domain_sizes_Mbp <- rbind(domain_sizes_Mbp, c("danRer11", "cluster", abs(zf_5DOM_3 - zf_3DOM_5)))

domain_order <- c("assembly", "5DOM", "cluster", "3DOM")
assembly_order <- c("mm39", "danRer11")

domain_sizes_Mbp$domain <- factor(domain_sizes_Mbp$domain, levels = domain_order)
domain_sizes_Mbp$assembly <- factor(domain_sizes_Mbp$assembly, levels = assembly_order)

# Plot mm39/danRer11

domain_sizes_Mbp <- domain_sizes_Mbp %>%
  mutate(size = as.numeric(size))

domain_sizes_Mbp_toPlt <- domain_sizes_Mbp %>%
  group_by(domain) %>%
  summarise(ratio = size[assembly == "mm39"] / size[assembly == "danRer11"])

y_intercept_value <- domain_sizes_Mbp_toPlt$ratio[domain_sizes_Mbp_toPlt$domain == "assembly"]

size_ratios <- ggplot(domain_sizes_Mbp_toPlt, aes(x = domain, y = ratio, color = domain)) +
  geom_point(size = 1.8) +
  geom_hline(yintercept = y_intercept_value, color = "red", linewidth = 0.1) +
  scale_color_manual("domain", values = c("#A9A9A9", "#006699", "#BE7FE6", "#339933")) +
  xlab("") +
  ylab("mm39 / danRer11") +
  theme_bw() +
  theme(legend.position = c(0.85, 0.85)) +
  theme(legend.position = "none", element_line(linewidth = 0.1), axis.text.x = element_text(angle = 45, hjust = 1, size = 7), axis.text.y = element_text(size = 7), axis.title = element_text(size = 7))

ggsave(paste0(size_ratios.pdfs, "1.pdf"), plot = size_ratios, units = "cm", width = 3.7, height = 5.5)

# Plot 5DOM/3DOM

domain_sizes_Mbp_toPlt2 <- domain_sizes_Mbp %>%
  group_by(assembly) %>%
  summarise(ratio = size[domain == "5DOM"] / size[domain == "3DOM"])

size_ratios_53 <- ggplot(domain_sizes_Mbp_toPlt2, aes(x = assembly, y = ratio)) +
  geom_point(size = 1.8) +
  geom_hline(yintercept = 1, color = "red", linewidth = 0.1) +
  xlab("") +
  ylab("5DOM / 3DOM") +
  theme_bw() +
  theme(legend.position = "none", element_line(linewidth = 0.1), axis.text.x = element_text(angle = 45, hjust = 1, size = 7), axis.text.y = element_text(size = 7), axis.title = element_text(size = 7))
size_ratios_53

ggsave(paste0(size_ratios.pdfs, "2.pdf"), plot = size_ratios_53, units = "cm", width = 2.4, height = 5.5)
