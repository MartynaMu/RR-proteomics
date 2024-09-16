# PCA -----------------------------------------------------------------
library(PCAtools)

ds = zscored
is.norm = "Medians-centered-cell-lines-zscored"
pca <- pca(drop_na(ds), metadata = annot_col, scale = TRUE)

p <- biplot(pca, 
       x = "PC1",
       y = "PC2",
       lab = NULL,
       shape = "Condition.L3",
       shapekey = c("CFPAC" = 15, "MiaPaca" = 17, "PANC1" = 8),
       showLoadings = FALSE,
       boxedLoadingsNames = FALSE, 
       colby = "Condition.L1",
       colkey = annot_colors$Condition.L1,
       title = is.norm,
       # encircle = TRUE,
       legendPosition = "right",
       hline = 0, 
       vline = 0)
filename = paste0("pca-", is.norm, ".png")
ggsave(p, filename=filename, path = "figures/allruns/final_quant/normalization/", width = 7, height = 5.5, dpi=100, scale=1, bg = "white")

screeplot <- screeplot(pca, title = is.norm, axisLabSize = 10)
filename = paste0("screeplot-", is.norm, ".png")
ggsave(screeplot, filename=filename, path = "figures/allruns/final_quant/normalization/", width = 7, height = 5, dpi=100, scale=1, bg = "white")

loadings <- plotloadings(pca, rangeRetain = 0.01, title = is.norm, labSize = 3)
filename = paste0("top1loadings-",is.norm,".png")
ggsave(loadings, filename=filename, path = "figures/allruns/final_quant/normalization/", width = 9, height = 12, dpi=100, scale=1, bg = "white")


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
