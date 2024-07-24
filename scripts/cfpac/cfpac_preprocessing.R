# 2D VS 3D YOUNG VS 3D OLD VS PDX - collaboration with RR
# Downstream processing of DIANN output
# INPUT: DIANN quant - FASTA with shared HUMAN and MOUSE proteins 
# CFPAC cell line
# Martyna Muszczek, 04.06.23

#Libraries
library(tidyverse)
library(diann)
library(smartPep)

#Import
pg_matrix <- read.delim("data/report.pg_matrix_allruns.tsv")

#store protein info in a separate df
protein_info <- pg_matrix[1:5]
df <- pg_matrix %>% select(matches("Genes|_CFPAC_"))

#tidy up colnames
temp <- merge2reps(df)

colnames(temp)[-1] <- colnames(df)[-1] %>%
  str_extract(pattern = "2D.+|3D.+|3Y.+|3O.+|Xeno.+|X.+") %>%
  str_replace_all(pattern = "_[:digit:].mzML.dia|_[:digit:]_repeat.+", replacement = "") %>%
  unique()

colnames(temp)[1] <- "Gene"

#long format
df_long <- pivot_longer(temp, cols=2:16, names_to = "Sample", values_to = "Intensity", values_drop_na = FALSE)

#Explores NAs
summary <- df_long %>% group_by(Sample) %>% summarise("NAs" = sum(is.na(Intensity)), 
                                                      "nonNAs" = sum(!is.na(Intensity)),
                                                      "%Missing" = NAs/nonNAs*100)
#correlation matrixes
corrplot::corrplot(cor(log2(drop_na(temp[-1]))),
                   method = "color",
                   hclust.method = "ward.D",
                   is.corr = FALSE)

#Heatmap global
pheatmap::pheatmap(drop_na(log2(temp[-1])), 
                   scale = "row", 
                   clustering_distance_cols = "euclidean",
                   show_rownames = FALSE)

boxplot(log2(df[-1])) # tech rep distribution
boxplot(log2(temp[-1])) # biol rep distribution

df_wide <- drop_na(log2(temp[-1]))

