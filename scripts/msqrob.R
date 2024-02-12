#Libs----------------------------------------------------------------------------
library("tidyverse")

# Load----------------------------------------------------------------------------
msqrob <- read.delim("data/msqrob_culture_vs_xeno.tsv")

#FC values
msqrob_diff <- msqrob[,19:26]
#Intensity values
msqrob <- msqrob[1:20]

# Tidying up--------------------------------------------------------------------
#drop nas
msqrob <- drop_na(msqrob)

# Tidy up gene groups - getting rid of mouse part and setting as rownames
msqrob$Gene <- str_replace(msqrob$Gene, pattern = "_HUMAN;.+|_HUMAN", replacement = "")
msqrob <- column_to_rownames(msqrob, "Gene")

# Correlation plots---------------------------------------------------------------
library(corrplot)

# AFTER MAXLFQ BEFORE MSQROB
corr <- cor(log2(maxlfq), use="complete.obs")
corrplot(corr, method = "color", hclust.method = "ward.D2", is.corr = FALSE, col = COL1("OrRd"))

# AFTER MSQROB
corr <- cor(msqrob[1:18])
corrplot(corr, method = "color", hclust.method = "ward.D2", is.corr = FALSE, col = COL1("OrRd"))

# AFTER QNORM
corr <- cor(qnorm[1:18])
corrplot(corr, method = "color", hclust.method = "complete", is.corr = FALSE, col = COL1("OrRd"))

# Global heatmap-----------------------------------------------------------------------
annot_col <- data.frame(
  Condition.L1 = as.factor(
    c(rep(c("2D", "3D.young", "3D.old"), each = 4),
      rep(c("PDX.2D", "PDX.3D"), each = 3))),
  Condition.L2 = as.factor(
    c(rep("2D", 4), 
      rep("3D", 8), 
      rep("PDX", 6))))

row.names(annot_col) <- colnames(msqrob)[1:18]

library(RColorBrewer)
annot_colors <- list(Condition.L1 = rev(brewer.pal(5, "Paired")),
                     Condition.L2 = c("firebrick1", "darkgreen", "blue"))

names(annot_colors$Condition.L1) <- levels(annot_col$Condition.L1)
names(annot_colors$Condition.L2) <- levels(annot_col$Condition.L2)

pheatmap::pheatmap(msqrob[1:18],
           scale = "row",
           show_rownames = FALSE,
           clustering_distance_cols = "euclidean",
           clustering_distance_rows = "correlation",
           annotation_col = annot_col,
           annotation_colors = annot_colors,
           main = "Global heatmap")

#as.dendrogram(hm$tree_col) %>% plot()

#display.brewer.pal(5, "Paired")
#brewer.pal(5, "Paired")

# PCA----------------------------------------------------------------------------
library(PCAtools)
#column to rownames for pca function
mat <- as.matrix(msqrob[1:18])
#pca
p <- pca(mat, metadata = annot_col, scale = TRUE)

png("figures/pca.png")
biplot(p, 
       lab = colnames(mat), 
       labSize = 5, 
       showLoadings = FALSE,
       boxedLoadingsNames = FALSE, 
       colby = "Condition.L1",
       colkey = annot_colors$Condition.L1,
       title = "Pre-normalisation",
       encircle = TRUE,
       hline = 0, 
       vline = 0)
dev.off()

screeplot(p)

#Norm methods-------------------------------------------------------------------
library(proBatch)
mat <- as.matrix(mat)
##median norm=========================
temp <- normalize_data_dm(mat, normalize_func ="medianCentering")
boxplot(temp)
##quantile norm=========================
qnorm <- normalize_data_dm(mat, normalize_func ="quantile")
rownames(qnorm) <- rownames(msqrob)
qnorm$Gene <- rownames(msqrob)
boxplot(temp)
##PARK7_HUMAN norm============
park7_mean <- mean(mat["PARK7_HUMAN",])
vec <- mat["PARK7_HUMAN",] - park7_mean
temp <- mat - vec
##scale normalization=============
temp <- scale(mat)
boxplot(temp)
##VSN norm=====================
temp <- limma::normalizeVSN(mat)
boxplot(temp)
##Loess norm=====================
temp <- limma::normalizeCyclicLoess(mat)
boxplot(temp)

##Assess PCA grouping---------------------------------------------------------
#pca
n <- pca(qnorm[1:18], metadata = annot_col, scale = TRUE)

png("figures/pca.png")
biplot(n, 
       lab = colnames(temp), 
       labSize = 5, 
       showLoadings = FALSE,
       boxedLoadingsNames = FALSE, 
       colby = "Condition.L1",
       colkey = annot_colors$Condition.L1,
       title = "Quantile normalisation",
       encircle = TRUE,
       hline = 0,
       vline = 0)
dev.off()

screeplot(n)

# Corr after norm--------------------------------------------------------------
# filter(msqrob, grepl("FN1", Gene))
qnorm <- as.data.frame(qnorm)
qnorm$Gene <- msqrob$Gene
rownames(qnorm) <- msqrob$Gene
pheatmap::pheatmap(qnorm[1:18],
                   scale = "row",
                   show_rownames = FALSE,
                   clustering_distance_cols = "correlation",
                   clustering_distance_rows = "correlation",
                   annotation_col = annot_col,
                   annotation_colors = annot_colors,
                   main = "Cell-cell signaling genes")

# GO classification------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
go_class <- groupGO(qnorm$Gene, 
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 level = 5, 
                 keyType = "SYMBOL")
go_class <- go_class@result

term = "cellular response to hypoxia"
filtr <- str_split_1(go_class$geneID[go_class$Description == term], pattern = "/")

# Corr of distinct GO terms --------------------------------------------------

pheatmap::pheatmap(qnorm[qnorm$Gene %in% filtr,1:18],
                   scale = "row",
                   show_rownames = TRUE,
                   clustering_distance_cols = "euclidean",
                   clustering_distance_rows = "correlation",
                   annotation_col = annot_col,
                   annotation_colors = annot_colors,
                   main = term)
