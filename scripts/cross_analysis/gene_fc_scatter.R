cell.line <- c("miapaca", "panc1", "cfpac")

logFCs <- data.frame(matrix(ncol = 1))
colnames(logFCs) <- "ID"

for (l in cell.line) {
  library(tidyverse)
  fc_msqrob <- read.delim(sprintf("data/%s/a_Hurdle_Msqrob.tsv",l))
  # fc_msqrob <- read.delim("data/panc1/a_msqrob_result_Ridge_Aggregate.tsv")
  
  # Preprocessing-------------------------------------------------------------------
  #tidy Gene symbols
  fc_msqrob$Genes <- str_replace_all(fc_msqrob$Genes, "_HUMAN", "")
  colnames(fc_msqrob) <- str_replace_all(colnames(fc_msqrob), "norm_", "")
  
  prot_info <- fc_msqrob %>% dplyr::select(matches("Gene|Protein|Peptide")) # protein info
  int_mat <- fc_msqrob %>% dplyr::select(!matches("FC|qval|_OR|_fisher|_pval")) # intensity matrix
  colnames(int_mat)
  
  #first check filtered on normal qval, then fisher_qval
  fc_msqrob <- fc_msqrob %>% dplyr::select(matches("FC|(?<!fisher)_qval|Genes", perl = TRUE)) 
  #fc_msqrob <- fc_msqrob %>% dplyr::select(matches("FC|(?<=fisher)_qval|Genes", perl = TRUE)) 
  fc_msqrob <- column_to_rownames(fc_msqrob, "Genes")
  colnames(fc_msqrob)
  
  fc_msqrob <- fc_msqrob[,order(colnames(fc_msqrob))]
  colnames(fc_msqrob)
  
  # get rid of 0
  int_mat[int_mat == 0] <- NA
  
  # normalize quantiles
  int_mat <- column_to_rownames(int_mat, "Genes")
  int_mat <- int_mat %>% dplyr::select(!matches("Protein|Peptide"))
  colnames(int_mat)
  qnorm <- limma::normalizeQuantiles(drop_na(int_mat))
  
  # DE limma --------------------------------------------------------------------
  # limma instead of hurdle/msqrob ttest
  
  library(tidyverse)
  library(limma)
  
  # In-depth comparisons-----------------------------------------------------------
  ##Design table creation-------------------------------------------------------------
  #where first digit in rep represents the column it will be filled with, second is nr of replicates
  
  mat <- as.matrix(qnorm[1:ncol(qnorm)])
  V1 <- str_count(colnames(mat), pattern = "2D(?!_XENO)") %>% sum()
  V2 <- str_count(colnames(mat), pattern = "3D_YOUNG") %>% sum()
  V3 <- str_count(colnames(mat), pattern = "3D_OLD") %>% sum()
  V4 <- str_count(colnames(mat), pattern = "2D_XENO") %>% sum()
  V5 <- str_count(colnames(mat), pattern = "3D_XENO") %>% sum()
  
  design <- model.matrix(~ 0+factor(c(rep(1,V1),
                                      rep(2,V2),
                                      rep(3,V3),
                                      rep(4,V4),
                                      rep(5,V5))))
  colnames(design) <- c("x2D", "x3D_young", "x3D_old", "PDX_2D", "PDX_3D")
  row.names(design) <- colnames(mat)
  
  ##Limma fit and matrix contrasts-------------------------------------------------
  fit <- lmFit(mat, design)
  cont.matrix <- makeContrasts(x2D_x3Dyoung = x3D_young - x2D, #1
                               x2D_x3Dold = x3D_old - x2D, #2
                               x2D_PDX2D = PDX_2D - x2D, #3
                               x3Dyoung_x3Dold = x3D_old - x3D_young, #4
                               x3Dold_PDX3D = PDX_3D - x3D_old, #5
                               x3Dyoung_xPDX3D = PDX_3D - x3D_young, #6
                               levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  
  fit_bayes <- eBayes(fit2) 
  
  for (i in coefs) {
    logFCs <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(int_mat)) %>%
              select(c(ID, logFC, adj.P.Val)) |>
              full_join(logFCs, by = "ID", suffix = names(coefs[i]))
  }
}