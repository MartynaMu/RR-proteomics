# CFPAC cell line experiment only
#Libs----------------------------------------------------------------------------
library("tidyverse")

#Normalize----------------------------------------------------------------------
##quantile norm=========================
qnorm <- limma::normalizeQuantiles(mat)

# Correlation plots---------------------------------------------------------------
library(corrplot)

corr <- cor(df_wide)
corrplot(corr, method = "color", hclust.method = "complete", is.corr = FALSE, col = COL1("OrRd"))

# AFTER QNORM
corr <- cor(qnorm)
corrplot(corr, method = "color", hclust.method = "complete", is.corr = FALSE, col = COL1("OrRd"))

# Global heatmap-----------------------------------------------------------------------
cfpac.l1 <- rep(c("2D", "3D.old", "3D.young", "PDX.2D", "PDX.3D"), each = 3)
miapaca.l1 <- rep(c("2D", "3D.old", "3D.young", "PDX.2D", "PDX.3D"), each = 3)
panc1.l1 <- c(rep(c("2D", "3D.old", "3D.young"), each = 4), rep(c("PDX.2D", "PDX.3D"), each = 3))

cfpac.l2 <- c(rep("2D", 3), rep("3D", 6), rep("PDX", 6))
miapaca.l2 <- c(rep("2D", 3), rep("3D", 6), rep("PDX", 6))
panc1.l2 <- c(rep("2D", 4), rep("3D", 8), rep("PDX", 6))

celllines <- c(rep("CFPAC", 15), rep("MiaPaca", 15), rep("PANC1", 18))


annot_col <- data.frame(
  Condition.L1 = as.factor(
    c(cfpac.l1,miapaca.l1,panc1.l1)),
  Condition.L2 = as.factor(
    c(cfpac.l2,miapaca.l2,panc1.l2)),
  Condition.L3 = as.factor(celllines))

row.names(annot_col) <- colnames(df_wide)

library(RColorBrewer)
annot_colors <- list(Condition.L1 = rev(brewer.pal(5, "Paired")),
                     Condition.L2 = c("firebrick1", "darkgreen", "blue"),
                     Condition.L3 = c("black", "darkgrey", "lightgray"))

names(annot_colors$Condition.L1) <- levels(annot_col$Condition.L1)
names(annot_colors$Condition.L2) <- levels(annot_col$Condition.L2)
names(annot_colors$Condition.L3) <- levels(annot_col$Condition.L3)

pheatmap::pheatmap(df_wide,
                   scale = "row",
                   show_rownames = FALSE,
                   clustering_distance_cols = "correlation",
                   clustering_distance_rows = "euclidean",
                   annotation_col = annot_col,
                   annotation_colors = annot_colors,
                   main = "Global heatmap")

# PCA----------------------------------------------------------------------------
library(PCAtools)
#pca
p <- pca(mat, metadata = annot_col, scale = TRUE)

png("figures/allruns/pca_prenorm.png")
biplot(p, 
      lab = NULL, 
      # labSize = 5, 
      shape = "Condition.L3",
      shapekey = c("CFPAC" = 15, "MiaPaca" = 17, "PANC1" = 8),
       showLoadings = FALSE,
       boxedLoadingsNames = FALSE, 
       colby = "Condition.L1",
       colkey = annot_colors$Condition.L1,
       title = "Pre-normalisation",
      # encircle = TRUE,
      legendPosition = "right",
       hline = 0, 
       vline = 0)
dev.off()

screeplot(p)

##Assess PCA grouping---------------------------------------------------------
#pca
n <- pca(qnorm, metadata = annot_col, scale = TRUE)

png("figures/allruns/pca_postnorm.png")
biplot(n, 
       lab = NULL, 
       # labSize = 5, 
       shape = "Condition.L3",
       shapekey = c("CFPAC" = 15, "MiaPaca" = 17, "PANC1" = 8),
       showLoadings = FALSE,
       boxedLoadingsNames = FALSE, 
       colby = "Condition.L1",
       colkey = annot_colors$Condition.L1,
       title = "Post-normalisation",
       # encircle = TRUE,
       legendPosition = "right",
       hline = 0, 
       vline = 0)
dev.off()

screeplot(n)

# Corr after norm--------------------------------------------------------------
# filter(msqrob, grepl("FN1", Genes))
qnorm <- as.data.frame(qnorm)
qnorm$Genes <- rownames(df_wide)
rownames(qnorm) <- qnorm$Genes

pheatmap::pheatmap(qnorm,
                   scale = "row",
                   show_rownames = FALSE,
                   clustering_distance_cols = "correlation",
                   clustering_distance_rows = "euclidean",
                   annotation_col = annot_col,
                   annotation_colors = annot_colors,
                   title = "Global heatmap after normalisation")

# GO classification------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
go_class <- groupGO(rownames(qnorm), 
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    level = 5, 
                    keyType = "SYMBOL")
go_class <- go_class@result %>% filter(Count > 0)

term = "positive regulation of extracellular matrix organization"
filtr <- str_split_1(go_class$geneID[go_class$Description == term], pattern = "/")

# Corr of distinct GO terms --------------------------------------------------

pheatmap::pheatmap(qnorm[rownames(qnorm) %in% filtr,],
                   scale = "row",
                   show_rownames = TRUE,
                   clustering_distance_cols = "euclidean",
                   clustering_distance_rows = "euclidean",
                   annotation_col = annot_col,
                   annotation_colors = annot_colors,
                   main = term)

