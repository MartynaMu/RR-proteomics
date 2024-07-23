# Input: a_msqrob_result_Ridge_Aggregate.tsv - FC calculated by msqrob with 
# ridge regression (you may wish to evaluate this) and msqrobAggregate 
# (very sensitive for small changes of very well quantified abundant proteins)

library(tidyverse)
fc_msqrob <- read.delim("data/panc1/a_Hurdle_Msqrob.tsv")
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

# IF USING MSQROB/HURDLE FC: change the sign to '-' so the reference is on the left
# fc_msqrob <- mutate_at(fc_msqrob, .vars = vars(ends_with("_FC")), .funs = list({function(x) x*(-1)}))

# EDA---------------------------------------------------------------------------
p <- pivot_longer(int_mat, 1:18, names_to = "Sample", values_to = "Intensity") %>%
        ggplot(aes(x=Sample, y=Intensity))+
        geom_violin()+
        geom_boxplot(width=.3)+
        ggpubr::rotate_x_text(angle=45)
ggsave(plot = p, filename="boxplots.png", device = "png",path = "figures/panc1", dpi = 100, units = "px", width=800)


p <- pivot_longer(qnorm, 1:18, names_to = "Sample", values_to = "Intensity") %>%
  ggplot(aes(x=Sample, y=Intensity))+
  geom_violin()+
  geom_boxplot(width=.3)+
  ggpubr::rotate_x_text(angle=45)
ggsave(plot = p, filename="boxplots-post-norm.png", device = "png",path = "figures/panc1", dpi = 100, units = "px", width=800)

