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
  # Retrieve overlaps and convert to descriptions ------------------------------
  # Assign whether terms are upreg or downreg
  # create a named list of GO IDs
  library(GO.db)
  full.overlap.ids <- c()
  for (i in curr.comp) {
    temp <- full.overlap[grepl(paste0(i,".upreg"), names(full.overlap))]
    temp <- base::Reduce(intersect, temp) %>% list() |> set_names(paste0(i,".upreg"))
    full.overlap.ids <- c(full.overlap.ids, temp)
    temp <- full.overlap[grepl(paste0(i,".downreg"), names(full.overlap))]
    temp <- base::Reduce(intersect, temp) %>% list() |> set_names(paste0(i,".downreg"))
    full.overlap.ids <- c(full.overlap.ids, temp)
  }
  # Convert GO IDs to term descriptions
  full.overlap.terms <- sapply(full.overlap.ids, Term)

return(full.overlap.terms)
}

full.overlap <- term.overlap(coefs = coefs, cell.line = cell.line, ontology = "WP")




