#find DEPs sign between comparisons -------------------------------------------
anova <- topTable(fit_bayes, adjust = "BH", genelist = rownames(mat), resort.by = "adj.P.val", number = nrow(qnorm))

anova_filtr <- anova$ProbeID[anova$adj.P.Val <= 0.05]
temp <- filter(qnorm, !Gene %in% anova_filtr)

pheatmap::pheatmap(temp[1:18],
                   scale = "row",
                   show_rownames = FALSE,
                   cutree_rows = 5,
                   clustering_distance_cols = "correlation",
                   clustering_distance_rows = "correlation",
                   annotation_col = annot_col,
                   annotation_colors = annot_colors,
                   main = "ANOVA filtered")

##Assess PCA grouping---------------------------------------------------------
#pca
library(PCAtools)
a <- pca(temp[1:18], metadata = annot_col, scale = TRUE)

png("figures/pca_anova_filtered.png")
biplot(a,
       lab = colnames(temp[1:18]), 
       labSize = 5, 
       showLoadings = FALSE,
       boxedLoadingsNames = FALSE, 
       colby = "Condition.L1",
       colkey = annot_colors$Condition.L1,
       title = "ANOVA filtered",
       encircle = TRUE,
       hline = 0,
       vline = 0,
       colLegendTitle = "Condition",
       legendPosition = "right",
       legendLabSize = 16)
dev.off()

screeplot(a)

# plot just legend from pca
library(ggpubr)
get_legend(temp) %>% as_ggplot()
