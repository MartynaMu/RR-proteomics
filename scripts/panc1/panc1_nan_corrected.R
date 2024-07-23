# Input: a_msqrob_result_Ridge_Aggregate.tsv - FC calculated by msqrob with 
# ridge regression (you may wish to evaluate this) and msqrobAggregate 
# (very sensitive for small changes of very well quantified abundant proteins)

library(tidyverse)
nan_corr <- read.delim("data/panc1/empty_corrected_for_visualizations_all.tsv")

#tidy Gene symbols
nan_corr$Genes <- str_replace_all(nan_corr$Genes, "_HUMAN", "")
colnames(nan_corr) <- str_replace_all(colnames(nan_corr), "norm_", "")

prot_info <- nan_corr %>% select(matches("Gene|Protein|Peptide")) # protein info
int_mat <- nan_corr %>% select(!matches("FC|qval|_OR|_fisher|_pval")) # intensity matrix
colnames(int_mat)

#first check filtered on normal qval, then fisher_qval
nan_corr <- nan_corr %>% select(matches("FC|(?<!fisher)_qval|Genes", perl = TRUE)) 
#nan_corr <- nan_corr %>% select(matches("FC|(?<=fisher)_qval|Genes", perl = TRUE)) 
nan_corr <- column_to_rownames(nan_corr, "Genes")
colnames(nan_corr)

nan_corr <- nan_corr[,order(colnames(nan_corr))]
colnames(nan_corr)
