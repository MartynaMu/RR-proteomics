# CFPAC cell line experiment only
#Libs----------------------------------------------------------------------------
library("tidyverse")

# Tidying up--------------------------------------------------------------------
#drop nas
mat <- as.matrix(df_wide)
##quantile norm=========================
qnorm <- limma::normalizeQuantiles(mat)

# Tidy up gene groups - getting rid of mouse part and setting as rownames
df_wide$Gene <- str_replace_all(df_wide$Gene, pattern = "_HUMAN", replacement = "")
df_wide <- column_to_rownames(df_wide, "Gene")

# Correlation plots---------------------------------------------------------------
library(corrplot)

corr <- cor(df_wide[1:18])
corrplot(corr, method = "color", hclust.method = "ward.D2", is.corr = FALSE, col = COL1("OrRd"))

# AFTER QNORM
corr <- cor(qnorm[1:18])
corrplot(corr, method = "color", hclust.method = "complete", is.corr = FALSE, col = COL1("OrRd"))

# Global heatmap-----------------------------------------------------------------------
annot_col <- data.frame(
  Condition.L1 = as.factor(
    c(rep(c("2D", "3D.old", "3D.young", "PDX.2D", "PDX.3D"), each = 3))),
  Condition.L2 = as.factor(
    c(rep("2D", 3), 
      rep("3D", 6), 
      rep("PDX", 6))))

row.names(annot_col) <- colnames(df_wide)

library(RColorBrewer)
annot_colors <- list(Condition.L1 = rev(brewer.pal(5, "Paired")),
                     Condition.L2 = c("firebrick1", "darkgreen", "blue"))

names(annot_colors$Condition.L1) <- levels(annot_col$Condition.L1)
names(annot_colors$Condition.L2) <- levels(annot_col$Condition.L2)

pheatmap::pheatmap(df_wide,
                   scale = "row",
                   show_rownames = FALSE,
                   clustering_distance_cols = "euclidean",
                   clustering_distance_rows = "correlation",
                   annotation_col = annot_col,
                   annotation_colors = annot_colors,
                   main = "Global heatmap")

# PCA----------------------------------------------------------------------------
library(PCAtools)
#column to rownames for pca function

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

#Normalization-------------------------------------------------------------------
##quantile norm=========================
qnorm <- limma::normalizeQuantiles(mat)

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
qnorm$Gene <- rownames(df_wide)
rownames(qnorm) <- qnorm$Gene

pheatmap::pheatmap(qnorm,
                   scale = "row",
                   show_rownames = FALSE,
                   clustering_distance_cols = "euclidean",
                   clustering_distance_rows = "correlation",
                   annotation_col = annot_col,
                   annotation_colors = annot_colors)

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
