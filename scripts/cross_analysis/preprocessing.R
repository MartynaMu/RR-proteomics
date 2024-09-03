library(tidyverse)
df <- read.delim("data/all_lines/a_Hurdle_Msqrob.tsv", header = TRUE, sep = "\t")
df <- df[1:49]

# Tidy gene names-----------------------------------------------
df$Genes <- str_replace_all(df$Genes, "_HUMAN", replacement = "")
df <- column_to_rownames(df, var = "Genes")

# Normalization ---------------------------------------------------
#qnorm <- limma::normalizeQuantiles(as.matrix(df))
qnorm <- limma::normalizeMedianValues(as.matrix(df))

# Subtract means-----------------------------------------------------
# identify columns of each cell line
panc.ind <- colnames(qnorm) %>% str_which(pattern="PANC")
miapaca.ind <- colnames(qnorm) %>% str_which(pattern="MIAPACA")
cfpac.ind <- colnames(qnorm) %>% str_which(pattern="CFPAC")

# calculate means of each cell line
mat <- qnorm %>% 
  as.data.frame() %>% 
  rowwise() %>% 
  mutate(PANC.mean = mean(c_across(all_of(panc.ind)),na.rm=TRUE),
         MIAPACA.mean = mean(c_across(all_of(miapaca.ind)), na.rm = TRUE),
         CFPAC.mean = mean(c_across(all_of(cfpac.ind)), na.rm = TRUE))

# subtract the means of cell lines from each
mat <- mutate(mat, across(.cols = all_of(panc.ind), .fns = function(x) (x - PANC.mean)),
               across(.cols = all_of(miapaca.ind), .fns = function(x) (x - MIAPACA.mean)),
               across(.cols = all_of(cfpac.ind), .fns = function(x) (x - CFPAC.mean)))
mat <- mat %>% as.data.frame() %>% drop_na()
mat <- mat[1:48]
rownames(mat) <- df %>% drop_na() %>% rownames()

# HM annotations ------------------------------------------
# Complex hm with annotations
cfpac.l1 <- rep(c("2D", "3D.young", "3D.old", "PDX.2D", "PDX.3D"), each = 3)
miapaca.l1 <- rep(c("2D", "3D.young", "3D.old", "PDX.2D", "PDX.3D"), each = 3)
panc1.l1 <- c(rep(c("2D", "3D.young", "3D.old"), each = 4), rep(c("PDX.2D", "PDX.3D"), each = 3))

cfpac.l2 <- c(rep("2D", 3), rep("3D", 6), rep("PDX", 6))
miapaca.l2 <- c(rep("2D", 3), rep("3D", 6), rep("PDX", 6))
panc1.l2 <- c(rep("2D", 4), rep("3D", 8), rep("PDX", 6))

celllines <- c(rep("PANC1", 18), rep("MiaPaca", 15), rep("CFPAC", 15))


annot_col <- data.frame(
  Condition.L1 = as.factor(
    c(panc1.l1, miapaca.l1, cfpac.l1)),
  Condition.L2 = as.factor(
    c(panc1.l2, miapaca.l2, cfpac.l2)),
  Condition.L3 = as.factor(celllines))

row.names(annot_col) <- colnames(mat[1:48])

library(RColorBrewer)
annot_colors <- list(Condition.L1 = rev(brewer.pal(5, "Paired")),
                     Condition.L2 = c("firebrick1", "darkgreen", "blue"),
                     Condition.L3 = c("black", "darkgrey", "lightgray"))

names(annot_colors$Condition.L1) <- levels(annot_col$Condition.L1)
names(annot_colors$Condition.L2) <- levels(annot_col$Condition.L2)
names(annot_colors$Condition.L3) <- levels(annot_col$Condition.L3)

# HM -----------------------------------------------------------------
pheatmap::pheatmap(mat,
                   annotation_col = annot_col,
                   annotation_colors = annot_colors,
                   annotation_names_col = FALSE,
                   show_rownames = FALSE,
                   main = "Global heatmap, medians normalized, means subtracted")

# PCA -----------------------------------------------------------------
library(PCAtools)
p <- pca(drop_na(df), metadata = annot_col, scale = TRUE) # pre norm
p <- pca(mat, metadata = annot_col, scale = TRUE) # post norm

png("figures/allruns/final_quant/pca-post-norm.png")
biplot(p, 
       x = "PC1",
       y = "PC2",
       lab = NULL, 
       # labSize = 5, 
       shape = "Condition.L3",
       shapekey = c("CFPAC" = 15, "MiaPaca" = 17, "PANC1" = 8),
       showLoadings = FALSE,
       boxedLoadingsNames = FALSE, 
       colby = "Condition.L1",
       colkey = annot_colors$Condition.L1,
       title = "Post-norm",
       # encircle = TRUE,
       legendPosition = "right",
       hline = 0, 
       vline = 0)
dev.off()

screeplot(p)

plotloadings(p)

