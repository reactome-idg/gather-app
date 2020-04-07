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
meta.data.path <- "/opt/data/GTEx/source/meta-data/"
expression.path <- "/opt/data/GTEx/source/raw-data/"

sample.attribute <- read.delim(paste0(meta.data.path, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"))
rnaseq.samples <- sample.attribute[which(sample.attribute$SMAFRZE %in% 'RNASEQ'), ]
rnaseq.samples <- rnaseq.samples[which(rnaseq.samples$SMRIN > 6.0), ]

tissue.summary <- rnaseq.samples %>% group_by(SMTSD) %>% tally(sort = T)
tissue.summary  <- tissue.summary[which(tissue.summary$n >= 30), ] # greater than 30 samples

rnaseq.samples <- rnaseq.samples[which(rnaseq.samples$SMTSD %in% tissue.summary$SMTSD), ]

gene_counts <- read_delim(paste0(expression.path, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"), delim = '\t', skip = 2)
gene_counts <- gene_counts[, c(1, 2, which(colnames(gene_counts) %in% rnaseq.samples$SAMPID))]
gene_counts$Name <- gsub("\\..*", "", gene_counts$Name)

# -------------------------------------------------------
# save all data as tissue specific matrix 56200 X 15185 
#-----------------------------------------------------
out.path.raw <- "/opt/data/GTEx/processed/raw-count/"
tissue.summary$file_name <- gsub('[)]', '', gsub('[(]', '-', gsub(' ', '', tissue.summary$SMTSD)))

lapply(tissue.summary$SMTSD, function(x){
  ids <- rnaseq.samples$SAMPID[which(rnaseq.samples$SMTSD %in% x)]
  dat <- as.data.frame(gene_counts[ , c(1, 2, which(colnames(gene_counts) %in% ids))])
  
  name <- tissue.summary$file_name[which(tissue.summary$SMTSD %in% x)]
  write.csv(dat, paste0(out.path.raw, paste0(name, ".csv")), row.names = T, quote = F)
})

#-----------------------------------------------------
# remove genes/engs that don't map to a unique location on gene
# CPM transform 
#-----------------------------------------------------
gene_counts <- gene_counts[which(!gene_counts$Description %in% unique(gene_counts$Description[duplicated(gene_counts$Description)])), ]
out.path.cpm <- "/opt/data/GTEx/processed/cpm/cpm/"
out.path.eda <- "/opt/data/GTEx/processed/cpm/eda/cpm-transformed/"
  
lapply(tissue.summary$SMTSD, function(tx){
  ids <- rnaseq.samples$SAMPID[which(rnaseq.samples$SMTSD %in% tx)]
  df <- as.data.frame(gene_counts[ , c(1, 2, which(colnames(gene_counts) %in% ids))])
  rownames(df) <- gene_counts$Name
  df <- df[, -c(1, 2)]
  
  tissue <- tissue.summary$file_name[which(tissue.summary$SMTSD %in% tx)]

  # samples are columns 
  y <- DGEList(count=df)
  x <- calcNormFactors(y, method="TMM")
  df.cpm <- cpm(x, normalized.lib.sizes=F, log=F)
  df.cpm <- as.data.frame(df.cpm)
  
  pdf(file = paste0(out.path.eda, paste0(tissue, ".pdf")),  wi = 12, he = 9)
  boxplot(log2(df.cpm[c(2000:3000), ]), las = 2, cex = 0.7, main = "", col = "lightblue", ylab = "log2(cpm(tmm))") 
  boxplot(log2(df[c(2000:3000), ]), las = 2, cex = 0.7, main = "", col = "lightblue", ylap = "log2(raw count)") 
  dev.off()
  
  file.name <- paste0(tissue, "_CPM.csv")
  write.csv(df.cpm, paste0(out.path.cpm, file.name), row.names = T, quote = F)
})

# ---------------------------------------------
# filter by gk_central 
# ----------------------------------------------
gene.subset <- gene_counts$Description[unique(gene_counts$Description) %in% gk_central$genes]
gene_counts <- gene_counts[which(gene_counts$Description %in% gene.subset), ]        
ensg_gene <- gene_counts[, c(1, 2)] 
names(ensg_gene) <- c("ensg", "gene")

write_csv(ensg_gene, "/opt/data/GTEx/processed/cpm/ensg_gene_all_final.csv")

# -----------------------------------------------
# make fi relation for correlation 
# -----------------------------------------------
fi.subset <- fi %>% filter(fi$from %in% gene_counts$Description & fi$to %in% gene_counts$Description)
fi.subset.genes <- unique(union(fi.subset$from, fi.subset$to))

g.subset <- gene_counts[, c(1, 2)]
colnames(g.subset) <- c('ensg_from', 'from')
df1 <- inner_join(fi.subset, g.subset, by = 'from')
colnames(g.subset) <- c('ensg_to', 'to')
df2 <- inner_join(df1,  g.subset, by = 'to')

df2 <- df2[!duplicated(df2[,1:2]), ]
length(union(unique(df2$from), unique(df2$to))) # 7691 total genes 
final.fi <- df2

write_csv(final.fi, "/opt/data/GTEx/processed/fi/ensg_gene_fi.csv")

# ------------------------------------------------
# do EDA over all genes in final adjacency 
# ------------------------------------------------
in.path <- "/opt/data/GTEx/processed/cpm/cpm/"
files <- list.files(path = in.path, pattern = "*.csv")

out.path.wo <- '/opt/data/GTEx/processed/cpm/wo-outliers/'
out.path.eda <- '/opt/data/GTEx/processed/cpm/eda/cpm/'


lapply(seq_along(files), function(i){
  f <- files[i]
  df <- as.data.frame(read_csv(paste0(in.path, f)))
  rownames(df) <- df[,1]
  df <- as.data.frame(df[ ,-1])
  df <- df[which(rownames(df) %in% ensg_gene$ensg), ]
  
  tissue <- gsub("_CPM.csv", "", f)
  
  df.meta <- rnaseq.samples[which(rnaseq.samples$SAMPID %in% names(df)), ]
  df.meta <- df.meta[match(df.meta$SAMPID, colnames(df)), ]
  
  df1 <- df
  df2 <- log2(df)
  df <- log2(df + 1)
  # sum(rowSums(df) == 0)
  
  PC <- prcomp(df, scale.=T, center = T)
  # Plot first 2 PCs
  plotdata <- data.frame(SampleID=rownames(PC$rotation),
                         PC1=PC$rotation[,1],
                         PC2=PC$rotation[,2])
  
  plotdata$batchid <- as.factor(df.meta$SMGEBTCH)
  
  samplepca.outlier <- tag.pca.outlier(plotdata$PC1, sampleqc = rownames(PC$rotation))
  
  g1 <- ggplot(plotdata, aes(x=PC1, y=PC2, group=batchid)) +
    geom_point(aes(color=rownames(PC$rotation),
                   text = paste("SampleID: ", rownames(PC$rotation),
                                "</br> outlier_tag: ", samplepca.outlier$pca_outlier))) +
    labs(title = tissue)
  p1 <- plotly::ggplotly(g1)
  
  write.csv(samplepca.outlier, paste0(out.path.eda, paste0(tissue, "-sample-outliers.csv")), row.names = F)
  htmlwidgets::saveWidget(as_widget(p1), paste0(out.path.eda, paste0(tissue, "-pca.html")), selfcontained = TRUE)
  
  df <- df[ ,which(!colnames(df) %in% samplepca.outlier$SampleID[which(samplepca.outlier$pca_outlier %in% '1')])]
  write.csv(df, paste0(out.path.wo, paste0(tissue, "_CPM.csv")), row.names = T, quote = F)
  write.csv(samplepca.outlier, paste0(out.path.eda, paste0(tissue, "-sample-outliers.csv")), row.names = F)
  
  pdf(file = paste0(out.path.eda, paste0(tissue, ".pdf")),  wi = 12, he = 9)
  boxplot(df1, las = 2, cex = 0.7, main = tissue, col = "lightblue", ylab = "cpm(tmm)")
  boxplot(df2, las = 2, cex = 0.7, main = tissue, col = "lightblue", ylab = "log2(cpm(tmm))")
  boxplot(df, las = 2, cex = 0.7, main = tissue, col = "lightblue", ylab = "log2(cpm(tmm) + 1) ")
  pheatmap::pheatmap(df[c(1000:2000), ])
  dev.off()
  
})

# ---------------------------------------
# RIN summary of tissues processed 
# ----------------------------------------
ggplot(rnaseq.samples, aes(x = reorder(SMTSD, -SMRIN, FUN = median), y = SMRIN)) + 
  geom_boxplot(fill = "lightblue", colour = "darkblue", outlier.alpha = 0.1) + 
  geom_hline(yintercept=6, linetype="dashed",  color = "darkorange", size=0.5) + 
  theme(text = element_text(size=15), axis.text.x = element_text(angle = 75, hjust = 1)) + 
  xlab('') + 
  ylab('RIN') + 
  labs(title = 'Summary of RIN RNA-Seq groupedby tissue collection site', subtitle = "Sample count > 30 & RIN > 6")
ggsave("/opt/data/GTEx/processed/cpm/RIN.png", width = 20, height = 15)


