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
ds = zscored
is.norm = "Medians-centered-cell-lines-zscored"
scale = "none"

hm <- pheatmap::pheatmap(ds,
                         scale = scale,
                         annotation_col = annot_col,
                         annotation_colors = annot_colors,
                         annotation_names_col = FALSE,
                         show_rownames = FALSE,
                         main = paste0("Global heatmap - ", is.norm),
                         col=colorRampPalette(c("blue","white","red"))(100))
filename = paste0("hm-",is.norm,".png")
ggsave(hm, filename=filename, path="figures/allruns/final_quant/normalization/", device="png", width = 8, height = 8, dpi=100, scale=1, bg = "white")

