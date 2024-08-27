library(tidyverse)
# Cross analysis of deps in all cell lines
coefs <- seq.int(1,18,1)
names(coefs) <- colnames(fit_bayes$contrasts)

cell.line <- c("PANC", "MIAPACA", "CFPAC")
  
# DF Creation------------------------------------------------------------------
# Read the deps files, merge per cell line, add a column with a comparison
for (l in cell.line) {
  files <- list.files(paste0("data/all_lines/comparison_stats/"), full.names = TRUE)
  temp <- files[str_which(files, pattern = l)] |> lapply(read_tsv)
  curr.coef.names <- names(coefs)[str_which(names(coefs), pattern = l)]
  names(temp) <- str_sort(curr.coef.names)
  temp <- mapply(c, temp, str_sort(curr.coef.names), SIMPLIFY = FALSE)
  temp <- data.table::rbindlist(temp)
  
  colnames(temp)[ncol(temp)] <- "Comparison"
 # temp <- temp %>% filter(!between(logFC,-1.3,1.3), adj.P.Val <= 0.05)
  temp <- temp %>% filter(adj.P.Val <= 0.05)
  assign(paste("deps",l, sep = "."), temp)
}

# Upregulated ------------------------------------------------------------------
## Create a list of elements to overlap with ------------------------------------
upreg <- c()
for (l in cell.line) {
  deps <- paste0("deps.", l)
  curr.coefs <- coefs[str_which(names(coefs), pattern = l)]
  for (i in curr.coefs) {
   # temp <- get(e) %>% filter(Category == "Up-regulated", Comparison == eval(names(coefs[i]))) %>% select(ID)
    temp <- get(deps) %>% filter(logFC > 0, Comparison == eval(names(curr.coefs[curr.coefs==i]))) %>% dplyr::select(ID)
    colnames(temp) <- paste("deps", l ,names(curr.coefs[curr.coefs==i]),sep=".")
    upreg <- c(upreg, temp)
  }
}

## Overlap ----------------------------------------------------------------------
library(eulerr)
curr.comp <- str_replace_all(names(coefs), pattern = "^[:alpha:]+.", replacement = "") %>% unique()
for (i in curr.comp) {
  comb <- upreg[grepl(i, names(upreg))]
  temp <- plot(euler(comb), quantities = TRUE, 
               labels = c("PANC", "MIAPACA", "CFPAC"), 
               main = i,
               fills = c("white", "lightgrey", "#C4E311"))
  ggsave(
    filename = paste0("upreg.deps.", i, ".png"),
    plot = temp,
    path = "figures/allruns/final_quant/overlap",
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
    # temp <- get(e) %>% filter(Category == "Up-regulated", Comparison == eval(names(coefs[i]))) %>% select(ID)
    temp <- get(deps) %>% filter(logFC < 0, Comparison == eval(names(curr.coefs[curr.coefs==i]))) %>% dplyr::select(ID)
    colnames(temp) <- paste("deps", l ,names(curr.coefs[curr.coefs==i]),sep=".")
    downreg <- c(downreg, temp)
  }
}

## Overlap ----------------------------------------------------------------------
library(eulerr)
curr.comp <- str_replace_all(names(coefs), pattern = "^[:alpha:]+.", replacement = "") %>% unique()
for (i in curr.comp) {
  comb <- downreg[grepl(i, names(downreg))]
  temp <- plot(euler(comb), quantities = TRUE, 
               labels = c("PANC", "MIAPACA", "CFPAC"), 
               main = i,
               fills = c("white", "lightgrey", "#C4E311"))
  ggsave(
    filename = paste0("downreg.deps.", i, ".png"),
    plot = temp,
    path = "figures/allruns/final_quant/overlap",
    device = "png",
    units = "px", 
    dpi = 100,
    width = 500,
    height = 500
  )
}

# Table ------------------------------------------------------------------------
library(gtExtras)
library(fontawesome)
temp <- bind_rows(deps.CFPAC,deps.MIAPACA,deps.PANC)
temp$Category <- NA
temp$Category[temp$logFC < 0] <- "Down-regulated"
temp$Category[temp$logFC > 0] <- "Up-regulated"
temp <- separate(temp, Comparison, into = c("Cell.line", "Comparison"), sep = "\\.")

temp <- temp %>% group_by(Cell.line, Comparison, Category) %>% summarise(n=n()) |> pivot_wider(names_from = "Comparison", values_from = "n")
temp <- mutate(temp, logo=ifelse(Category == "Up-regulated", 
                                 "arrow-up", 
                                 "arrow-down"))

temp %>% 
  gt() %>% 
  gt_theme_espn() |>
  cols_merge(c(Category, logo)) |>
  fmt_icon(logo,fill_color = c("arrow-up" = "red", "arrow-down" = "blue")) %>%
  cols_label(logo="")
