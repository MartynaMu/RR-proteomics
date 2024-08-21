library(tidyverse)
# Cross analysis of deps in all cell lines
coefs <- seq.int(1,18,1)
names(coefs) <- colnames(fit_bayes$contrasts)

cell.line <- c("PANC", "MIAPACA", "CFPAC")

term.overlap <- function(coefs, cell.line, ontology = "BP|CC|MF|KEGG|WP") {
  # DF Creation------------------------------------------------------------------
  # Read the deps files, merge per cell line, add a column with a comparison
  if (ontology %in% c("BP", "CC", "MF")) {
    for (l in cell.line) {
      files <- list.files(paste0("figures/allruns/final_quant/gsego_results/gsego_", ontology, "_results/"), 
                          full.names = TRUE,
                          pattern = "unsimplified")
      
      temp <- files[str_which(files, pattern = l)] |> lapply(read_tsv)
      curr.coef.names <- names(coefs)[str_which(names(coefs), pattern = l)]
      names(temp) <- str_sort(curr.coef.names)
      temp <- mapply(c, temp, str_sort(curr.coef.names), SIMPLIFY = FALSE)
      temp <- data.table::rbindlist(temp)
      
      colnames(temp)[ncol(temp)] <- "Comparison"
      assign(paste("deps",l, sep = "."), temp)
    }
  } else if (ontology == "KEGG") {
    for (l in cell.line) {
      files <- list.files(paste0("figures/allruns/final_quant/gsekegg_results/"), 
                          full.names = TRUE,
                          pattern = "unsimplified")
      
      temp <- files[str_which(files, pattern = l)] |> lapply(read_tsv)
      curr.coef.names <- names(coefs)[str_which(names(coefs), pattern = l)]
      names(temp) <- str_sort(curr.coef.names)
      temp <- mapply(c, temp, str_sort(curr.coef.names), SIMPLIFY = FALSE)
      temp <- data.table::rbindlist(temp)
      
      colnames(temp)[ncol(temp)] <- "Comparison"
      assign(paste("deps",l, sep = "."), temp)
    } 
  } else if (ontology == "WP") {
    for (l in cell.line) {
      files <- list.files(paste0("figures/allruns/final_quant/gseWP_results/"), 
                          full.names = TRUE,
                          pattern = ".tsv")
      
      temp <- files[str_which(files, pattern = l)] |> lapply(read_tsv)
      curr.coef.names <- names(coefs)[str_which(names(coefs), pattern = l)]
      names(temp) <- str_sort(curr.coef.names)
      temp <- mapply(c, temp, str_sort(curr.coef.names), SIMPLIFY = FALSE)
      temp <- data.table::rbindlist(temp)
      
      colnames(temp)[ncol(temp)] <- "Comparison"
      assign(paste("deps",l, sep = "."), temp)
    }
  }
  # Upregulated ------------------------------------------------------------------
  ## Create a list of elements to overlap with ------------------------------------
  full.overlap <- c()
  
  upreg <- c()
  for (l in cell.line) {
    deps <- paste0("deps.", l)
    curr.coefs <- coefs[str_which(names(coefs), pattern = l)]
    for (i in curr.coefs) {
      temp <- get(deps) %>% filter(NES > 0, Comparison == eval(names(curr.coefs[curr.coefs==i]))) %>% dplyr::select(ID)
      colnames(temp) <- paste("deps", l ,names(curr.coefs[curr.coefs==i]),sep=".")
      upreg <- c(upreg, temp)
    }
  }
  names(upreg) <- str_c(names(upreg), ".upreg")
  full.overlap <- c(full.overlap, upreg)
  
  ## Overlap ----------------------------------------------------------------------
  library(eulerr)
  curr.comp <- str_replace_all(names(coefs), pattern = "^[:alpha:]+.", replacement = "") %>% unique()
  for (i in curr.comp) {
    comb <- upreg[grepl(i, names(upreg))]
    temp <- plot(euler(comb), 
                 quantities = TRUE, 
                 labels = c("PANC", "MIAPACA", "CFPAC"), 
                 main = i,
                 fills = c("white", "lightgrey", "#C4E311"))
    ggsave(
      filename = paste0("upreg.terms.", i, ".png"),
      plot = temp,
      path = paste0("figures/allruns/final_quant/overlap/terms",ontology),
      device = "png",
      units = "px", 
      dpi = 100,
      width = 500,
      height = 500
    )
  }
  
  # Down-regulated ------------------------------------------------------------------
  ## Create a list of elements to overlap with ------------------------------------
  downreg <- c()
  for (l in cell.line) {
    deps <- paste0("deps.", l)
    curr.coefs <- coefs[str_which(names(coefs), pattern = l)]
    for (i in curr.coefs) {
      temp <- get(deps) %>% filter(NES < 0, Comparison == eval(names(curr.coefs[curr.coefs==i]))) %>% dplyr::select(ID)
      colnames(temp) <- paste("deps", l ,names(curr.coefs[curr.coefs==i]),sep=".")
      downreg <- c(downreg, temp)
    }
  }
  names(downreg) <- str_c(names(downreg), ".downreg")
  full.overlap <- c(full.overlap, downreg)
  
  ## Overlap ----------------------------------------------------------------------
  library(eulerr)
  curr.comp <- str_replace_all(names(coefs), pattern = "^[:alpha:]+.", replacement = "") %>% unique()
  for (i in curr.comp) {
    comb <- downreg[grepl(i, names(downreg))]
    temp <- plot(euler(comb), 
                 quantities = TRUE, 
                 labels = c("PANC", "MIAPACA", "CFPAC"),
                 main = i,
                 fills = c("white", "lightgrey", "#C4E311"))
    ggsave(
      filename = paste0("downreg.terms.", i, ".png"),
      plot = temp,
      path = paste0("figures/allruns/final_quant/overlap/terms",ontology),
      device = "png",
      units = "px", 
      dpi = 100,
      width = 500,
      height = 500
    )
  }
  return(full.overlap)
}

full.overlap <- term.overlap(coefs = coefs, cell.line = cell.line, ontology = "WP")

# Retrieve overlaps of GO and convert to descriptions ------------------------------
# Assign whether terms are upreg or downreg
# doesn't work if there is no full overlap (check diagrams)
# create a named list of GO IDs
# library(GO.db)
# full.overlap.ids <- c()
# for (i in curr.comp) {
#   temp <- full.overlap[grepl(paste0(i,".upreg"), names(full.overlap))]
#   temp <- purrr::reduce(temp, intersect) %>% list() |> set_names(paste0(i,".upreg"))
#   full.overlap.ids <- c(full.overlap.ids, temp)
#   temp <- full.overlap[grepl(paste0(i,".downreg"), names(full.overlap))]
#   temp <- purrr::reduce(temp, intersect) %>% list() |> set_names(paste0(i,".downreg"))
#   full.overlap.ids <- c(full.overlap.ids, temp)
# }

# Convert GO IDs to term descriptions
#full.overlap.terms <- sapply(full.overlap.ids, Term)

# Term overlaps in at least 2 cell lines
full.overlap.up <- c()
for (i in curr.comp) {
  temp <- c()
  temp1 <- full.overlap[grepl(paste0(i,".upreg"), names(full.overlap))]
  temp <- c(temp,intersect(temp1[[1]], temp1[[2]]))
  temp <- c(temp,intersect(temp1[[1]], temp1[[3]]))
  temp <- c(temp,intersect(temp1[[3]], temp1[[2]]))
  temp <- list(temp)
  full.overlap.up <- c(full.overlap.up,temp)
}

names(full.overlap.up) <- paste0(curr.comp,".upreg")

full.overlap.down <- c()
for (i in curr.comp) {
  temp <- c()
  temp1 <- full.overlap[grepl(paste0(i,".down"), names(full.overlap))]
  temp <- c(temp,intersect(temp1[[1]], temp1[[2]]))
  temp <- c(temp,intersect(temp1[[1]], temp1[[3]]))
  temp <- c(temp,intersect(temp1[[3]], temp1[[2]]))
  temp <- list(temp)
  full.overlap.down <- c(full.overlap.down,temp)
}

names(full.overlap.down) <- paste0(curr.comp,".down")

full.overlap.terms <- c(full.overlap.down,full.overlap.up)

# Retrieve gene symbols from overlapping terms
# 1. prepare a list of pathways enriched in at least 2 cell lines
temp <- unlist(full.overlap.terms) %>% unique()

# 2. Find each pathway in gsea objects and retrieve genes from the pathway
term.genes <- c()
for (i in wp_list) {
  term.genes <- c(term.genes, geneInCategory(i)[i$ID %in% temp])
}

term.desc <- c()
for (i in wp_list) {
  term.desc <- c(term.desc, i[i$ID %in% temp]$Description)
}
term.desc <- term.desc %>% unique()

wp_desc <- data.frame("ID" = temp, "Description" = term.desc)
# 3. store as a vector to use for search in matrix
term.genes.vec <- unlist(term.genes) %>% unique() %>% sort()

# 4. Prepare a matrix of FC and p-values from each comparison in each cell line
args <- c(number = nrow(mat), resort.by = "logFC")

temp <- data.frame("Genes" = rownames(topTable(fit_bayes, number=nrow(mat))))
for (i in seq.int(1,6,1)) {
  temp <- left_join(temp,(topTable(fit_bayes,coef=i,substitute(args))) |> select(logFC) %>% rownames_to_column(var = "Genes"), by="Genes")
  temp <- left_join(temp,(topTable(fit_bayes,coef=i+6,substitute(args))) |> select(logFC) %>% rownames_to_column(var = "Genes"), by="Genes")
  temp <- left_join(temp,(topTable(fit_bayes,coef=i+12,substitute(args))) |> select(logFC) %>% rownames_to_column(var = "Genes"), by="Genes")
  colnames(temp)[(ncol(temp)-2):ncol(temp)] <- str_c(cell.line,curr.comp[i],"logFC", sep = ".")
}
assign("FC", temp)
FC <- column_to_rownames(FC, var="Genes")

# 5. Retrieve FC of the genes from each comparison in each cell line
FC_wp <- FC %>% filter(rownames(FC) %in% term.genes.vec)

# 6. Look for a specific term genes
term <- term.genes[["WP2363"]]


# 7. Draw heatmap
p <- pheatmap::pheatmap(FC_wp[term,],
                   cluster_cols = FALSE,
                   gaps_col = c(3,6,9,12,15),
                   cluster_rows = TRUE,
                   angle_col = 45,
                   display_numbers = TRUE,
                   fontsize_number = 10,
                   number_color = "black",
                   border_color = "black",
                   main = "Log2FC in each comparison",
                   col=colorRampPalette(c("cornflowerblue","white","red"))(100))

ggsave(filename = "hm_wp_pathway.png",
       plot = p,
      device = "png",
      height = 500,
      width = 1100,
      units = "px",
      dpi = 100)

# Loop for all terms --------------------------------------------------------------
terms <- names(term.genes)
for (i in terms) {
  if (length(term.genes[[i]]) > 20) {
    p <- pheatmap::pheatmap(FC_wp[term.genes[[i]],], 
                            cluster_cols = FALSE,
                            gaps_col = c(3,6,9,12,15),
                            cluster_rows = TRUE,
                            angle_col = 45,
                            border_color = "black",
                            main = "Log2FC in each comparison",
                            col=colorRampPalette(c("cornflowerblue","white","red"))(100))
    ggsave(filename = paste0(i,".hm.png"),
           plot = p,
           path = "figures/allruns/final_quant/term_overlap/",
           device = "png",
           height = 500,
           width = 1100,
           units = "px",
           dpi = 100)
  } else if (length(term.genes[[i]]) <= 20) {
    p <- pheatmap::pheatmap(FC_wp[term.genes[[i]],], 
                            cluster_cols = FALSE,
                            gaps_col = c(3,6,9,12,15),
                            cluster_rows = TRUE,
                            angle_col = 45,
                            display_numbers = TRUE,
                            fontsize_number = 10,
                            number_color = "black",
                            border_color = "black",
                            main = "Log2FC in each comparison",
                            col=colorRampPalette(c("cornflowerblue","white","red"))(100))
    ggsave(filename = paste0(i,".hm.png"),
           plot = p,
           path = "figures/allruns/final_quant/term_overlap/",
           device = "png",
           height = 500,
           width = 1100,
           units = "px",
           dpi = 100)
  }
}


