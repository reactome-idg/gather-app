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
p_load("edgeR")
# ----------------------------------------------
args = commandArgs(trailingOnly=T)
# print(args)

if(length(args) == 0) {
   print("Need to pass paths to input/output dirs.")
}
consortia <- args[1]

if (consortia = "tcga") {
    in.path <- "/opt/data/TCGA/processed/cpm/wo-outliers/"
    out.path <- "/opt/data/TCGA/processed/fi/cpm-spearman/"
    fi.path <- "/opt/data/TCGA/processed/fi/ensg_gene_fi.csv"

    in.path.wel <- "/opt/data/TCGA/processed/fi/cpm-spearman/"
    out.path.wel <- "/opt/data/TCGA/processed/fi/wel-spearman/"
}

if (consortia = "gtex") {
    in.path <- "/opt/data/GTEx/processed/cpm/wo-outliers/"
    out.path <- "/opt/data/GTEx/processed/fi/cpm-spearman/"
    fi.path <- "/opt/data/GTEx/processed/fi/ensg_gene_fi.csv"

    in.path.wel <- "/opt/data/GTEx/processed/fi/cpm-spearman/"
    out.path.wel <- "/opt/data/GTEx/processed/fi/wel-spearman/"
}
# -----------------------------------------------
options(stringsAsFactors = FALSE)
set.seed(1234)
source("functions.R", encoding = "UTF-8")
# ----------------------------------------------------------------
# make Functional Interaction (FI) spearman rho rank
# ----------------------------------------------------------------
dir.ls <- list.files(path = in.path, pattern = ".csv")

lapply(seq_along(dir.ls), function(num){
  f <- dir.ls[num]
  tissue <- gsub("_CPM.csv", "", f)

  df <- as.data.frame(read_csv(paste0(in.path, f)))
  rownames(df) <- df[, 1]
  df <- df[, -1]

  fi <- read_csv("/opt/data/TCGA/processed/fi/ensg_gene_fi.csv")

  sp <- lapply(seq_along(1:dim(fi)[1]), function(i){
    cor(as.numeric(df[which(rownames(df) %in% fi$ensg_from[i]), ]), as.numeric(df[which(rownames(df) %in% fi$ensg_to[i]), ]), method ="spearman")
  })
  sp <- unlist(sp)
  fi$weight <-  sp

  write_csv(fi,  paste0(paste0(out.path, tissue), "_Spearman_FI.csv"), col_names = T)
})

# --------------------------------------------------------------
# weighted edge list file  formats for mcl input
# --------------------------------------------------------------
dir.ls <- list.files(path = in.path.wel, pattern = "*.csv")

lapply(seq_along(dir.ls), function(i){
  f <- dir.ls[i]
  tissue <- gsub("_Spearman_FI.csv", "", f)

  dat <- read_csv(paste0(in.path.wel, f))
  dat <- dat[, c("from", "to", "weight")]
  dat <- dat[complete.cases(dat$weight), ]

  dat$weight <- abs(dat$weight)

  write_delim(dat, paste0(out.path.wel, paste0(tissue, "_Spearman.txt")), delim = '\t', col_names = F)
})
