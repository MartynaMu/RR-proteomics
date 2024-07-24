library(tidyverse)
# Cross analysis of deps in all cell lines
coefs <- seq.int(1,6,1)
names(coefs) <- c("2Dv3Dy", "2Dv3Do", "2DvPDX2D", "3Dyv3Do", "3DovPDX3D", "3DyvPDX3D")

cell.line <- c("panc1", "cfpac", "miapaca")

# DF Creation------------------------------------------------------------------
# Read the deps files, merge per cell line, add a column with a comparison
for (l in cell.line) {
  files <- list.files(paste0("figures/", l, "/int_mat_fc/gsego_results/gsego_BP_results/"), 
                      full.names = TRUE,
                      pattern = ".tsv")
  temp <- lapply(files, read_tsv)
  names(temp) <- str_sort(names(coefs))
  temp <- mapply(c, temp, str_sort(names(coefs)), SIMPLIFY = FALSE)
  temp <- data.table::rbindlist(temp)
  
  colnames(temp)[ncol(temp)] <- "Comparison"
  assign(paste("deps",l, sep = "."), temp)
}

# Upregulated ------------------------------------------------------------------
## Create a list of elements to overlap with ------------------------------------
deps <- lapply("deps", paste, cell.line, sep = ".")
upreg <- c()
for (i in coefs) {
  for (e in deps) {
    temp <- get(e) %>% filter(NES > 0, Comparison == eval(names(coefs[i]))) %>% select(ID)
    colnames(temp) <- paste(e,names(coefs[i]),sep=".")
    upreg <- c(upreg, temp)
  }
}

## Overlap ----------------------------------------------------------------------
library(eulerr)
for (i in coefs) {
  comb <- upreg[grepl(names(coefs[i]), names(upreg))]
  temp <- plot(euler(comb), quantities = TRUE, labels = c("CFPAC", "MiaPaca", "PANC1"), main = names(coefs[i]),
               fills = c("white", "lightgrey", "#C4E311"))
  ggsave(
    filename = paste0("upreg.terms.", names(coefs[i]), ".png"),
    plot = temp,
    path = "figures/cross-analysis/overlap/bp",
    device = "png",
    units = "px", 
    dpi = 100,
    width = 500,
    height = 500
  )
}

# Down-regulated ------------------------------------------------------------------
## Create a list of elements to overlap with ------------------------------------
deps <- c("deps.cfpac", "deps.miapaca", "deps.panc1")
downreg <- c()
for (i in coefs) {
  for (e in deps) {
    temp <- get(e) %>% filter(NES < 0, Comparison == eval(names(coefs[i]))) %>% select(ID)
    colnames(temp) <- paste(e,names(coefs[i]),sep=".")
    downreg <- c(downreg, temp)
  }
}

## Overlap ----------------------------------------------------------------------
library(eulerr)
for (i in coefs) {
  comb <- downreg[grepl(names(coefs[i]), names(downreg))]
  temp <- plot(euler(comb), quantities = TRUE, labels = c("CFPAC", "MiaPaca", "PANC1"), main = names(coefs[i]),
               fills = c("white", "lightgrey", "#C4E311"))
  ggsave(
    filename = paste0("downreg.terms.", names(coefs[i]), ".png"),
    plot = temp,
    path = "figures/cross-analysis/overlap/bp/",
    device = "png",
    units = "px", 
    dpi = 100,
    width = 500,
    height = 500
  )
}


# DEBUG ------------------------------------------------------------------------
