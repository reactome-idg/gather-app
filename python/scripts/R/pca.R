# -------------------------------------------------------------------------
# 3D PCA plots of GTEx and TCGA FI CPM expressions 
# -------------------------------------------------------------------------
usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
usePackage("pacman")
p_load("tidyverse")
p_load("foreach")
p_load("reshape")
p_load("dplyr")
p_load("plotly")
p_load("plyr")
p_load("edgeR")

# -------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
set.seed(1234)
# -------------------------------------------------------------------------
# GTEx
# -------------------------------------------------------------------------
idg.fi.path <- "/opt/output/idg_fi.csv"
meta.data.path <- "/opt/GTEx/source/meta-data/"
cpmwooutlier.path <- "/opt/GTEx/processed/cpm/wo-outliers/"
sample.attribute <- read.delim(paste0(meta.data.path, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"))
sample.attribute$tissue <-  gsub('[)]', '', gsub('[(]', '-', gsub(' ', '', sample.attribute$SMTSD)))

idg.fi <- read_csv(idg.fi.path)
files <- list.files(cpmwooutlier.path)

samples.wooutlier <- lapply(seq_along(files), function(i){
  f <- files[i]
  tissue <- gsub("_CPM.csv", "", f)
  
  df <- as.data.frame(read_csv(paste0(cpmwooutlier.path, f)))
  rownames(df) <- df[,1]
  df <- df[, -1]
  
  return(tibble(samples = colnames(df), tissue = tissue))
})
all.samples <- do.call(rbind, samples.wooutlier)

sample.attribute.wo <-  sample.attribute[which(sample.attribute$SAMPID %in% all.samples$samples), ]

# -------------------------------------------------------------------------
in.path <- "/opt/GTEx/processed/cpm/wo-outliers/"
out.path <- "/opt/GTEx/processed/cpm/"
mappings <- read_csv("/opt/GTEx/processed/cpm/ensg_gene_all_final.csv")

df.ls <- lapply(seq_along(files), function(i){
  f <- files[i]
  tissue <- gsub("_CPM.csv", "", f)
  
  df <- as.data.frame(read_csv(paste0(cpmwooutlier.path, f)))
  rownames(df) <- df[,1]
  df <- df[, -1]
  
  return(df)
})
all.df <- do.call(cbind, df.ls)
all.df <- all.df[, match(sample.attribute.wo$SAMPID, colnames(all.df))]
rownames(all.df) <- mappings$gene

all.df.fi <- all.df[which(rownames(all.df) %in% idg.fi$genes), ]
# dim(all.df.fi)
# [1]  7197 14918
# -------------------------------------------------------------------------
PC <- prcomp(all.df.fi, scale.=T, center = T)

PC1 <- PC$rotation[,1]
PC2 <- PC$rotation[,2]
PC3 <- PC$rotation[,3]

pc.var <- round(((PC$sdev^2) / (sum(PC$sdev^2)))*100, 2)[1:10]
p <- plot_ly() %>%
  add_trace(x=PC1, y=PC2, z=PC3,  color = ~sample.attribute.wo$SMTS, colors = "Set2",
            type="scatter3d", mode="markers", opacity = 1, marker = list(size = 4)) %>%
  layout(
    title = "GTEx Tissues",
    scene = list(
      xaxis = list(title = paste(paste("PC1 ", pc.var[1]), "%")),
      yaxis = list(title = paste(paste("PC2 ", pc.var[2]), "%")),
      zaxis = list(title = paste(paste("PC3 ", pc.var[3]), "%"))
    ))


htmlwidgets::saveWidget(p, "/opt/output/GTEx_pca_3d_2020_fi.html", selfcontained = TRUE)
# -------------------------------------------------------------------------
# TCGA
# -------------------------------------------------------------------------
in.path <- "/opt/TCGA/processed/cpm/wo-outliers/"
count.barcode <-  read_csv("/opt/TCGA/source/meta-data/count_deduplicated_barcode.csv")
mappings <- read_csv("/opt/TCGA/processed/cpm/ensg_gene_all_final.csv")

dir.ls <- list.files(path = in.path, pattern = "*.csv")

df.ls <- lapply(seq_along(dir.ls), function(i){
  f <- dir.ls[i]
  tissue <- gsub("_CPM.csv", "", f)
  
  df <- as.data.frame(read_csv(paste0(in.path, f)))
  rownames(df) <- df[,1]
  df <- df[, -1]
  
  
  return(df)
})
all.df <- do.call(cbind, df.ls)

count.barcode <- count.barcode[which(count.barcode$submitter_id %in% names(all.df)), ]
count.barcode <- count.barcode[match(colnames(all.df), count.barcode$submitter_id), ]

rownames(all.df) <- mappings$gene
all.df.fi <- all.df[which(rownames(all.df) %in% idg.fi$genes), ]
# > dim(all.df.fi)
# [1] 7232 9495

PC <- prcomp(all.df.fi, scale.=T, center = T)

PC1 <- PC$rotation[,1]
PC2 <- PC$rotation[,2]
PC3 <- PC$rotation[,3]

pc.var <- round(((PC$sdev^2) / (sum(PC$sdev^2)))*100, 2)[1:10]
p <- plot_ly() %>%
  add_trace(x=PC1, y=PC2, z=PC3,  color = ~count.barcode$project.name, colors = "Set2",
            type="scatter3d", mode="markers", opacity = 1, marker = list(size = 4)) %>%
  layout(
    title = "TCGA Projects",
    scene = list(
      xaxis = list(title = paste(paste("PC1 ", pc.var[1]), "%")),
      yaxis = list(title = paste(paste("PC2 ", pc.var[2]), "%")),
      zaxis = list(title = paste(paste("PC3 ", pc.var[3]), "%"))
    )) 


htmlwidgets::saveWidget(p, "/opt/output/tcga_pca_3d_2020_fi.html", selfcontained = TRUE)