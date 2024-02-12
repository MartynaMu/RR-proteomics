#Libs----------------------------------------------------------------------------
library("tidyverse")

#Load------------
pg_matrix <- read.delim("data/report.pg_matrix.tsv")

#store protein info in a separate df
protein_info <- pg_matrix[1:5]
df <- pg_matrix[c(3, 6:41)]

df <- column_to_rownames(df, var = "Protein.Names")
temp <- colnames(df)
#merge
df <- smartPep::merge2reps(df)

colnames(df) <- str_sub(temp, start=38) %>% 
  str_replace_all(pattern = "_1.+$|_2.+$|_2_repeat.+$", replacement = "") %>% 
  unique() %>%
  str_c("x", .)

df <- log2(df)

#drop nas
df <- drop_na(df)

#pheatmap-----------------------------------------------------------------------
annot_col <- data.frame(
  Condition.L1 = as.factor(
    c(rep(c("2D", "3D.young", "3D.old"), each = 4),
      rep(c("PDX.2D", "PDX.3D"), each = 3))),
  Condition.L2 = as.factor(
    c(rep("2D", 4), 
      rep("3D", 8), 
      rep("PDX", 6))))

row.names(annot_col) <- colnames(df)[1:18]

library(RColorBrewer)
annot_colors <- list(Condition.L1 = rev(brewer.pal(5, "Paired")),
                     Condition.L2 = c("firebrick1", "darkgreen", "blue"))

names(annot_colors$Condition.L1) <- levels(annot_col$Condition.L1)
names(annot_colors$Condition.L2) <- levels(annot_col$Condition.L2)

hm <- pheatmap::pheatmap(df[1:18], 
                         scale = "row",
                         show_rownames = FALSE,
                         clustering_distance_cols = "correlation",
                         annotation_col = annot_col,
                         annotation_colors = annot_colors)

as.dendrogram(hm$tree_col) %>% plot()

#pca----------------------------------------------------------------------------
library(PCAtools)
#column to rownames for pca function
mat <- as.matrix(df)
#pca
p <- pca(mat, metadata = annot_col, scale = TRUE)

png("figures/pca.png")
biplot(p, 
       lab = colnames(df), 
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
##median norm=========================
temp <- normalize_data_dm(mat, normalize_func ="medianCentering")
boxplot(temp)
##quantile norm=========================
temp <- normalize_data_dm(mat, normalize_func ="quantile")
boxplot(temp)
##PARK7_HUMAN norm============
park7_mean <- mean(mat["PARK7_HUMAN",])
vec <- mat["PARK7_HUMAN",] - park7_mean
temp <- mat - vec
##scale normalization============
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
n <- pca(temp[,-c(5:8)], 
         metadata = annot_col, 
         scale = TRUE)

png("figures/pca.png")
biplot(n, 
       lab = colnames(temp[,-c(5:8)]), 
       labSize = 5, 
       showLoadings = FALSE,
       boxedLoadingsNames = FALSE, 
       colby = "Condition.L1",
       title = "Loess normalisation",
       colkey = annot_colors$Condition.L1,
       hline = 0,
       vline = 0,
       encircle = TRUE
)
dev.off()

screeplot(n)
