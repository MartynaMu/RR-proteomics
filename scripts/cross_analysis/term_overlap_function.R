library(tidyverse)
# Cross analysis of deps in all cell lines
coefs <- seq.int(1,6,1)
names(coefs) <- c("2Dv3Dy", "2Dv3Do", "2DvPDX2D", "3Dyv3Do", "3DovPDX3D", "3DyvPDX3D")

cell.line <- c("panc1", "cfpac", "miapaca")

term.overlap <- function(coefs, cell.line, ontology = "BP|CC|MF|KEGG") {
  # DF Creation------------------------------------------------------------------
  # Read the deps files, merge per cell line, add a column with a comparison
  if (ontology != "KEGG") {
    for (l in cell.line) {
      files <- list.files(paste0("figures/", l, "/int_mat_fc/gsego_results/gsego_", ontology, "_results/"), 
                          full.names = TRUE,
                          pattern = ".tsv")
      temp <- lapply(files, read_tsv, show_col_types=FALSE)
      names(temp) <- str_sort(names(coefs))
      temp <- mapply(c, temp, str_sort(names(coefs)), SIMPLIFY = FALSE)
      temp <- data.table::rbindlist(temp)
      
      colnames(temp)[ncol(temp)] <- "Comparison"
      assign(paste("deps",l, sep = "."), temp)
    }
  } else if (ontology == "KEGG") {
    for (l in cell.line) {
      files <- list.files(paste0("figures/", l, "/int_mat_fc/gsekegg_results/"), 
                          full.names = TRUE,
                          pattern = ".tsv")
      temp <- lapply(files, read_tsv, show_col_types=FALSE)
      names(temp) <- str_sort(names(coefs))
      temp <- mapply(c, temp, str_sort(names(coefs)), SIMPLIFY = FALSE)
      
      colnames(temp)[ncol(temp)] <- "Comparison"
      assign(paste("deps",l, sep = "."), temp)
    }
  }
  
  #create a list of variable names to access them later
  deps <- str_c("deps", cell.line, sep = ".")
  full.overlap <- c()
  
  # Upregulated ------------------------------------------------------------------
  ## Create a list of elements to overlap with ------------------------------------
  upreg <- c()
  for (i in coefs) {
    for (e in deps) {
      temp <- get(e) %>% filter(NES > 0, Comparison == eval(names(coefs[i]))) %>% dplyr::select(ID)
      colnames(temp) <- paste(e,names(coefs[i]),sep=".")
      upreg <- c(upreg, temp)
    }
  }
  names(upreg) <- str_c(names(upreg), ".upreg")
  full.overlap <- c(full.overlap, upreg)
  
  ## Overlap ----------------------------------------------------------------------
  library(eulerr)
  for (i in coefs) {
    comb <- upreg[grepl(names(coefs[i]), names(upreg))]
    temp <- plot(euler(comb), quantities = TRUE, labels = c("CFPAC", "MiaPaca", "PANC1"), main = names(coefs[i]),
                 fills = c("white", "lightgrey", "#C4E311"))
    ggsave(
      filename = paste0("upreg.terms.", names(coefs[i]), ".png"),
      plot = temp,
      path = paste0("figures/cross-analysis/overlap/",ontology),
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
  for (i in coefs) {
    for (e in deps) {
      temp <- get(e) %>% filter(NES < 0, Comparison == eval(names(coefs[i]))) %>% dplyr::select(ID)
      colnames(temp) <- paste(e,names(coefs[i]),sep=".")
      downreg <- c(downreg, temp)
    }
  }
  names(downreg) <- str_c(names(downreg), ".downreg")
  full.overlap <- c(full.overlap, downreg)
  
  ## Overlap ----------------------------------------------------------------------
  library(eulerr)
  for (i in coefs) {
    comb <- downreg[grepl(names(coefs[i]), names(downreg))]
    temp <- plot(euler(comb), quantities = TRUE, labels = c("CFPAC", "MiaPaca", "PANC1"), main = names(coefs[i]),
                 fills = c("white", "lightgrey", "#C4E311"))
    ggsave(
      filename = paste0("downreg.terms.", names(coefs[i]), ".png"),
      plot = temp,
      path = paste0("figures/cross-analysis/overlap/", ontology),
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
  for (i in coefs) {
    temp <- full.overlap[grepl(paste0(names(coefs[i]),".upreg"), names(full.overlap))]
    temp <- purrr::reduce(temp, intersect) %>% list() |> set_names(paste0(names(coefs[i]),".upreg"))
    full.overlap.ids <- c(full.overlap.ids, temp)
    temp <- full.overlap[grepl(paste0(names(coefs[i]),".downreg"), names(full.overlap))]
    temp <- purrr::reduce(temp, intersect) %>% list() |> set_names(paste0(names(coefs[i]),".downreg"))
    full.overlap.ids <- c(full.overlap.ids, temp)
  }
  # Convert GO IDs to term descriptions
  full.overlap.terms <- sapply(full.overlap.ids, Term)
  
return(full.overlap.terms)
}


full.overlap <- term.overlap(coefs = coefs, cell.line = cell.line, ontology = "BP")




