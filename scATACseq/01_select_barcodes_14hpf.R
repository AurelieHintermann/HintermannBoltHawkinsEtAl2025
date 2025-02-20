options(timeout = 10000)
# And package to read the metadata
if (!"openxlsx" %in% installed.packages()) {
install.packages("openxlsx", repos = "https://stat.ethz.ch/CRAN/")
}
library(openxlsx)  # Download the metadata

library(stringr) # Get last characters
target.directory <- "scATACseq/data/"
dir.create(target.directory, showWarnings = FALSE, recursive = TRUE)
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE243256&format=file&file=GSE243256%5FZEPA%5Fmetadata%2Exlsx"
if (!file.exists(file.path(target.directory, "GSE243256_ZEPA_metadata.xlsx"))) {
    download.file(url, file.path(target.directory, "GSE243256_ZEPA_metadata.xlsx"))
}
cell.metadata <- read.xlsx(xlsxFile = file.path(target.directory, "GSE243256_ZEPA_metadata.xlsx"))
cell.metadata$stage.integer <- as.numeric(gsub("hpf", "", cell.metadata$stage))
# cell.metadata$well <- str_sub(cell.metadata$cell_names, start= -4)
table(cell.metadata$stage.integer)
#      4      5      6      7      8      9     10     11     12     14     18 
#   9696  10203   9244  22487  23527  15751  17723  33078  16318  47518  73396 
#     20     22     24     30     34     38     42     48     72 
#  33079  72331 111471  60217  59929  43232  40175 101793  87645 
my.stage <- 14
cat(cell.metadata$cell_barcodes[cell.metadata$stage.integer == my.stage], sep = "\n", file = paste0("scATACseq/cell_barcodes_", my.stage, ".txt"))
print(table(cell.metadata$cell_type[cell.metadata$stage.integer == my.stage]))
#   hpf14:Anterior neural ridge            hpf14:Blood island 
#                          1863                           806 
#            hpf14:Diencephalon hpf14:Differentiating neurons 
#                           746                          1094 
#                hpf14:Doublets             hpf14:Endothelial 
#                          2314                           288 
#               hpf14:Epidermal     hpf14:Epidermal (foxi3a+) 
#                          3999                           223 
#              hpf14:Floorplate                     hpf14:Gut 
#                           327                           838 
#          hpf14:Hatching gland             hpf14:Heart field 
#                           746                          1444 
#        hpf14:Hindbrain dorsal       hpf14:Hindbrain ventral 
#                          1923                          2737 
# hpf14:Lateral line primordium  hpf14:Mesoderm_adaxial cells 
#                           851                           535 
#                hpf14:Midbrain        hpf14:Midbrain ventral 
#                          2920                           366 
#                 hpf14:Myotome            hpf14:Neural crest 
#                          2530                          1564 
#  hpf14:Neural plate posterior               hpf14:Notochord 
#                          2998                           922 
#        hpf14:Optic primordium            hpf14:Otic placode 
#                          2959                           469 
#                hpf14:Periderm         hpf14:Pharyngeal arch 
#                           459                           887 
#         hpf14:Pronephric duct        hpf14:Tailbud mesoderm 
#                           338                          3923 
#     hpf14:Tailbud spinal cord           hpf14:Telencephalon 
#                          3541                          2699 
#                     hpf14:YSL 
#                           209 
