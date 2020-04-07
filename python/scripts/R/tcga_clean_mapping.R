usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
usePackage("pacman")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

p_load("tidyverse")
p_load("TCGAutils")

# -----------------------------------------------
# duplication problem 
# https://www.biostars.org/p/311017/
# *file name and case id mappings
# -----------------------------------------------
count.md <- read_csv("/opt/data/TCGA/source/meta-data/tcga-api-count-metadata.csv")
count.entity.barcodes <- filenameToBarcode(count.md$file_name)
count.barcode <- left_join(count.md, count.entity.barcodes, by = c('file_id', 'file_name'))

colnames(count.barcode)[1] <- "submitter_id"
colnames(count.barcode)[55] <- "aliquot_submitter_id"

count.barcode <- count.barcode %>% 
  separate(aliquot_submitter_id, into = c('barcode_project', 'barcode_tss', 'barcode_participant' ,'barcode_sample_vial', 'barcode_portion_analyte', 'barcode_plate', 'barcode_center'), sep="-", remove = F)

count.barcode$short_barcode <- apply(count.barcode[, c('submitter_id', 'barcode_sample_vial')] , 1 , paste , collapse = "-" )
count.barcode <- count.barcode[order(count.barcode$case_id, count.barcode$barcode_plate, count.barcode$short_barcode), ]
count.barcode <- count.barcode[!duplicated(count.barcode$case_id), ]
write_csv(count.barcode, "/opt/data/TCGA/source/meta-data/count_deduplicated_barcode.csv")
