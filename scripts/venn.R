# Libraries---------------------------------------------------------------------
library(eulerr)
library(tidyverse)
library(limma)

coefs <- seq.int(1,6,1)
names(coefs) <- c("2Dv3Dy", "2Dv3Do", "2DvPDX2D", "3Dyv3Do", "3DovPDX3D", "3DyvPDX3D")


# Create euler object -------------------------------------------------------------
deps_venn <- function(fit){
  comb <- list()
  for (i in coefs) {
    deps <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(qnorm))
    
    volc_data <- deps
    volc_data$Category <- "Not significant"
    volc_data$Category[volc_data$logFC <= -1.3 & volc_data$adj.P.Val <= 0.05] <- "Down-regulated"
    volc_data$Category[volc_data$logFC >= 1.3 & volc_data$adj.P.Val <= 0.05] <- "Up-regulated"
    volc_data$Label <- NA
    volc_data$Label[volc_data$Category != "Not significant"] <- volc_data$ID[volc_data$Category != "Not significant"]
    
    upreg <- volc_data$Label[volc_data$Category == "Up-regulated"]
    downreg <- volc_data$Label[volc_data$Category == "Down-regulated"]

    comb <- append(comb, c(list(upreg), list(downreg)))
  }
  comb_names <- rep((names(coefs)),each=2,times=1) |> paste0(rep(c("_upreg", "_downreg"), every=2))
  names(comb) <- comb_names
  venn <- euler(comb)
  return(venn)
}

deps_venn(fit_bayes)

# Create list of DEPs of all groups ------------------------------------------
deps_list <- function(fit){
  comb <- list()
  for (i in coefs) {
    deps <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(qnorm))
    
    volc_data <- deps
    volc_data$Category <- "Not significant"
    volc_data$Category[volc_data$logFC <= -1.3 & volc_data$adj.P.Val <= 0.05] <- "Down-regulated"
    volc_data$Category[volc_data$logFC >= 1.3 & volc_data$adj.P.Val <= 0.05] <- "Up-regulated"
    volc_data$Label <- NA
    volc_data$Label[volc_data$Category != "Not significant"] <- volc_data$ID[volc_data$Category != "Not significant"]
    
    upreg <- volc_data$Label[volc_data$Category == "Up-regulated"]
    downreg <- volc_data$Label[volc_data$Category == "Down-regulated"]
    
    comb <- append(comb, c(list(upreg), list(downreg)))
  }
  comb_names <- rep((names(coefs)),each=2,times=1) |> paste0(rep(c("_upreg", "_downreg"), every=2))
  names(comb) <- comb_names
  return(comb)
}

deps_list <- deps_list(fit_bayes)

# Visualize on hm --------------------------------------------------------------
filtr_deps <- unlist(deps_list) |> unique()
temp <- pheatmap::pheatmap(qnorm[qnorm$Gene %in% filtr_deps,1:18],
                   scale = "row",
                   show_rownames = FALSE,
                   cutree_rows = 5,
                   annotation_row = clusters,
                   clustering_distance_cols = "euclidean",
                   clustering_distance_rows = "correlation",
                   annotation_col = annot_col,
                   annotation_colors = annot_colors,
                   main = "Significant DEPs")

# View row tree, decide on nr of clusters
temp$tree_row %>% 
  as.dendrogram() %>%
  plot(horiz = TRUE)

# assign genes to clusters and annotate in hm, extract and explore in string network/go
clusters <- cutree(temp$tree_row, k = 5)
clusters <- as.data.frame(clusters) |> mutate(clusters=as.factor(clusters))
clipr::write_clip(rownames(clusters)[clusters$clusters == 5])

