usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
usePackage("pacman")

p_load("tidyverse")
p_load("reshape")
p_load("dplyr")
p_load("plotly")
p_load("plyr")
p_load("pheatmap")
p_load("edgeR")

# -----------------------------------------------
options(stringsAsFactors = FALSE)
set.seed(1234)
source("functions.R", encoding = "UTF-8")
# -----------------------------------------------
fi <- read.delim("/opt/data/curated/FIs_Reactome_Genes.txt", header = F)
fi.genes <- unique(union(fi[,1], fi[,2]))
colnames(fi) <- c('from', 'to')

dgenes <- read.delim("/opt/data/curated/TCRD_Protein_Type_V_6_0_0.txt", header = T)

gk_central <- read_delim("/opt/data/curated/HumanGenesInReactome_091119.txt", delim = '\n' , col_names = F)
colnames(gk_central) <- c("genes")

# -----------------------------------------------
# total ensg genes in TCGA 60484
# -----------------------------------------------
in.path <- "/opt/data/TCGA/processed/raw-count/"
dir.ls <- list.files(path = in.path, pattern = "*.csv")

f <- dir.ls[1]
project <- strsplit(f, '_')[[1]][1]
gbm <- as.data.frame(read_csv(paste0(in.path, f)))
rownames(gbm) <- gbm[, 1]
gbm <- gbm[, -1]

# -----------------------------------------------
# 56553 out 60484 of mapped in ensemble - fetched from ui
# -----------------------------------------------
ens.mappping <- read_csv("/opt/data/TCGA/source/ui-fetch/ensemble-fetch/mart_export-unique-genes.txt")
names(ens.mappping)[1] <- "ensg"
names(ens.mappping)[5] <- "gene"

dat.ensg <- as.data.frame(rownames(gbm)[which(rownames(gbm) %in% unique(ens.mappping$ensg))])
colnames(dat.ensg) <- "ensg"
dim(dat.ensg)

dat.ensg.gene <- left_join(dat.ensg, ens.mappping[, c(1, 5)], by = 'ensg')
dup.genes <- unique(dat.ensg.gene$gene[duplicated(dat.ensg.gene$gene)])
dat.ensg.gene <- dat.ensg.gene[!duplicated(dat.ensg.gene), ]
dat.ensg.gene <- dat.ensg.gene[which(!dat.ensg.gene$gene %in% unique(dat.ensg.gene$gene[duplicated(dat.ensg.gene$gene)])), ]

# -----------------------------------------------
# transform count to cpm after deduplication 
# -----------------------------------------------
in.path <- "/opt/data/TCGA/processed/raw-count/"
out.path <- "/opt/data/TCGA/processed/cpm/cpm/"
out.path.eda  <- "/opt/data/TCGA/processed/cpm/eda/cpm-transformed/"

files <- list.files(path = in.path, pattern = "*.csv")

lapply(seq_along(files), function(i){
  df <- as.data.frame(read_csv(paste0(in.path, files[i])))
  rownames(df) <- df[ ,1]
  df <- df[,-1]
  df <- df[which(rownames(df) %in% dat.ensg.gene$ensg), ]
  df <- as.matrix(df)
  project <- strsplit(files[i], "_")[[1]][1]
  
  # samples are columns 
  y <- DGEList(count=df)
  x <- calcNormFactors(y, method="TMM")
  df.cpm <- cpm(x, normalized.lib.sizes=F, log=F)
  df.cpm <- as.data.frame(df.cpm)
  
  pdf(file = paste0(out.path.eda, paste0(project, ".pdf")),  wi = 12, he = 9)
  boxplot(log2(df.cpm[c(1:1000), ]), las = 2, cex = 0.7, main = "", col = "lightblue", ylab = "log2(cpm(tmm))") 
  boxplot(log2(df[c(1:1000), ]), las = 2, cex = 0.7, main = "", col = "lightblue", ylap = "log2(raw count)") 
  dev.off()
  
  file.name <- paste0(project, "_CPM.csv")
  write.csv(df.cpm, paste0(out.path, file.name), row.names = T, quote = F)
})

# ----------------------------------------------
# filter by gk_central 
# ----------------------------------------------
dat.ensg.gene <- dat.ensg.gene[which(dat.ensg.gene$gene %in% gk_central$genes), ]

write_csv(dat.ensg.gene, "/opt/data/TCGA/processed/cpm/ensg_gene_all_final.csv")

# -----------------------------------------------
# make fi relation for correlation 
# -----------------------------------------------
fi.subset <- fi %>% filter(fi$from %in% dat.ensg.gene$gene & fi$to %in% dat.ensg.gene$gene)
fi.subset.genes <- unique(union(fi.subset$from, fi.subset$to))

g.subset <- dat.ensg.gene
colnames(g.subset) <- c('ensg_from', 'from')
df1 <- inner_join(fi.subset, g.subset, by = 'from')
colnames(g.subset) <- c('ensg_to', 'to')
df2 <- inner_join(df1,  g.subset, by = 'to')

df2 <- df2[!duplicated(df2[,1:2]), ]
length(union(unique(df2$from), unique(df2$to))) # 7685 total genes 
final.fi <- df2
write_csv(final.fi, "/opt/data/TCGA/processed/fi/ensg_gene_fi.csv")

# ------------------------------------------------
# do EDA over all genes in final adjacency 
# ------------------------------------------------
in.path <- "/opt/data/TCGA/processed/cpm/cpm/"
out.path <- "/opt/data/TCGA/processed/cpm/eda/cpm/"
out.path.wo <- "/opt/data/TCGA/processed/cpm/wo-outliers/"
count.barcode <-  read_csv("/opt/data/TCGA/source/meta-data/count_deduplicated_barcode.csv")

dir.ls <- list.files(path = in.path, pattern = "*.csv")

lapply(seq_along(dir.ls), function(i){
  f <- dir.ls[i]
  project <- strsplit(f, '_')[[1]][1]
  gbm <- as.data.frame(read_csv(paste0(in.path, f)))
  rownames(gbm) <- gbm[, 1]
  gbm <- gbm[, -1]
  
  # filter by final genes to make entire adjacency total 19979 genes
  df <- gbm[which(rownames(gbm) %in% dat.ensg.gene$ensg), ]
  count.barcode <- count.barcode[match(colnames(df), count.barcode$submitter_id), ]
  
  df.trans1 <- df
  df.trans2 <- log2(df)
  df <- log2(df + 1)
  
  PC <- prcomp(df, scale.=T, center = T)
  # Plot first 2 PCs
  plotdata <- data.frame(SampleID=rownames(PC$rotation),
                         PC1=PC$rotation[,1],
                         PC2=PC$rotation[,2])
  
  plotdata$aliquot_submitter_id <- as.factor(count.barcode$aliquot_submitter_id)
  
  samplepca.outlier <- tag.pca.outlier(plotdata$PC1, sampleqc = rownames(PC$rotation))
  
  g1 <- ggplot(plotdata, aes(x=PC1, y=PC2, group=aliquot_submitter_id)) +
    geom_point(aes(color=aliquot_submitter_id,
                   text = paste("SampleID: ", rownames(PC$rotation),
                                "</br> outlier_tag: ", samplepca.outlier$pca_outlier))) +
    labs(title = project)
  p1 <- plotly::ggplotly(g1)
  write.csv(samplepca.outlier, paste0(out.path, paste0(project, "-sample-outliers.csv")), row.names = F)
  htmlwidgets::saveWidget(as_widget(p1), paste0(out.path, paste0(project, "-pca.html")), selfcontained = TRUE)
  
  df <- df[,which(!colnames(df) %in% samplepca.outlier$SampleID[which(samplepca.outlier$pca_outlier %in% '1')])]
  write.csv(df, paste0(out.path.wo, paste0(project, "_CPM.csv")), row.names = T, quote = F)
  write.csv(samplepca.outlier, paste0(out.path, paste0(project, "-sample-outliers.csv")), row.names = F)
  
  pdf(file = paste0(out.path, paste0(project, ".pdf")),  wi = 12, he = 9)
  boxplot(df.trans1, las = 2, cex = 0.7, main = project, col = "lightblue", ylab = "cpm(tmm)")
  boxplot(df.trans2, las = 2, cex = 0.7, main = project, col = "lightblue", ylab = "log2(cpm(tmm))")
  boxplot(df, las = 2, cex = 0.7, main = project, col = "lightblue", ylab = "log2(cpm(tmm) + 1) ")
  pheatmap::pheatmap(df[sample(nrow(df), 1000), ])
  dev.off()
})
