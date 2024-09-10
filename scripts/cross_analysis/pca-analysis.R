# PCA -----------------------------------------------------------------
library(PCAtools)
pca <- pca(drop_na(df), metadata = annot_col, scale = TRUE) # pre norm
pca <- pca(mat, metadata = annot_col, scale = TRUE) # post norm

biplot(pca, 
       x = "PC1",
       y = "PC2",
       lab = NULL,
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

screeplot(pca)

loadings <- plotloadings(pca, rangeRetain = 0.01,title = "Post-norm",labSize = 3)
ggsave(loadings, filename="figures/allruns/final_quant/loadingstop1.png", width = 8,height = 8, dpi=100, scale=1, bg = "white")


pairsplot <- pairsplot(pca,
          shape = "Condition.L3",
          shapekey = c("CFPAC" = 15, "MiaPaca" = 17, "PANC1" = 8),
          colby = "Condition.L1",
          colkey = annot_colors$Condition.L1,
          hline = 0, 
          vline = 0,
          triangle = TRUE, 
          trianglelabSize = 12,
          pointSize = 3,
          gridlines.major = FALSE, 
          gridlines.minor = FALSE,
          plotaxes = FALSE,
          title="Post-norm")
ggsave(pairsplot, filename="figures/allruns/final_quant/pairsplot.png", width = 8,height = 8, dpi=100, scale=1, bg = "white")
