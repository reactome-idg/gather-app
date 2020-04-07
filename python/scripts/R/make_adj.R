usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
usePackage("pacman")

p_load("tidyverse")
p_load("reshape")
p_load("plyr")

# -----------------------------------------------
options(stringsAsFactors = FALSE)
set.seed(1234)
source("functions.R", encoding = "UTF-8")
# ----------------------------------------------
args = commandArgs(trailingOnly=T)
# print(args)

if(length(args) == 0) {
   print("Need to pass paths to input/output dirs.")
}
consortia <- args[1]

if (consortia = "tcga") {
    in.path <- "/opt/data/TCGA/processed/cpm/wo-outliers/"
    out.path <- "/opt/data/TCGA/processed/adjacency/"
    mapping <- read_csv("/opt/data/TCGA/processed/cpm/ensg_gene_all_final.csv")
}

if (consortia = "gtex") {
    in.path <- "/opt/data/GTEx/processed/cpm/wo-outliers/"
    out.path <- "/opt/data/GTEx/processed/adjacency/"
    mapping <- read_csv("/opt/data/GTEx/processed/cpm/ensg_gene_all_final.csv")
}

# -----------------------------------------------
# make adjacency
# -----------------------------------------------
dir.ls <- list.files(in.path)
dir.ls.out <- list.files(out.path)

# if the run was't successful for all files redo 
if(length(list.files(path = out.path, pattern = "*.csv", full.names = TRUE)) != 0){
  in.tissue <- gsub("_CPM.csv", "", dir.ls)
  out.tissue <- gsub("_Spearman_Adj.csv", "", dir.ls.out)
  dir.ls <- dir.ls[which(!in.tissue %in% out.tissue)]
  
  # filter by date example
  # in.path <- "/Users/sanati/Box/reactome-nasim/gtex2/data/tissue/ea-spearman"
  # filter.date <- '2019-11-01'
  # 
  # details <- file.info(list.files(path = in.path, pattern = "*.csv", full.names = TRUE))
  # details <- details[with(details, order(as.POSIXct(mtime))), ]
  # details <- details[grepl(filter.date, as.POSIXct(details$mtime)), ]
  # files <- rownames(details)
}

lapply(seq_along(dir.ls), function(i){
  f <- dir.ls[i]
  tissue <- gsub("_CPM.csv", "", f)
  
  df <- as.data.frame(read_csv(paste0(in.path, f)))
  rownames(df) <- df[, 1]
  df <- df[, -1]
  
  df <- as.matrix(t(df))
  dat.cor <- rownames_to_column(as.data.frame(cor(df, method = "spearman")))

  # replace ensg with gene name wf
  # names(df) <- mapping$gene[match(names(df), mapping$ensg)]
  # rownames(df) <- mapping$gene[match(rownames(df), mapping$ensg)]
  # df <- as_tibble(as_tibble(rownames_to_column(df)))
  # names(df)[1] <- ""

  write_csv(dat.cor,  paste0(paste0(out.path, tissue), "_Spearman_Adj.csv"), col_names = T)
  rm(df, dat.cor, tissue) # could possible help with memory manag
  # check for logical values in NA rows
})

