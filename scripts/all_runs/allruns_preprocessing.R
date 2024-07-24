# 2D VS 3D YOUNG VS 3D OLD VS PDX - collaboration with RR
# Downstream processing of DIANN output
# INPUT: DIANN quant - FASTA with shared HUMAN and MOUSE proteins 
# Martyna Muszczek, 04.06.23

#Libraries
library(tidyverse)
library(diann)
library(smartPep)

#Import
# pg_matrix <- read.delim("data/report.pg_matrix_allruns.tsv")
df <- read.delim("data/report.gg_matrix_allruns.tsv")

#store protein info in a separate df
# protein_info <- pg_matrix[1:5]
# df <- pg_matrix[c(3, 6:103)]

# Removing outliers, merging -------------------------------------------------------
df <- df %>% 
  select(!matches("3O3_1|2D1_2_repeat"))

temp <- merge2reps(df)

# Tidying up rownames------------------------------------------------------------
temp$Genes <- str_replace_all(temp$Genes, pattern = "_HUMAN", replacement = "")
temp <- column_to_rownames(temp, "Genes")

# Tidying colnames----------------------------------------------------------------

colnames(temp) <- colnames(df)[-1] %>%
  str_extract(pattern = "2D.+|3D.+|3Y.+|3O.+|Xeno.+|X.+") %>%
  str_replace_all(pattern = "_[:digit:].mzML.dia|_[:digit:]_repeat.+", replacement = "") %>%
  unique()

colnames(temp)

colnames(temp)[str_detect(colnames(temp),"CFPAC")] <- str_c(rep("CFPAC_",15), 
                                                            c(rep("2D_",3), rep("3DO_", 3), 
                                                              rep("3DY_", 3), rep("PDX2D_",3), 
                                                              rep("PDX3D_", 3)), 
                                                            rep(c("1", "2", "3"),5))
colnames(temp)[str_detect(colnames(temp),"MiaPaca")] <- str_c(rep("MiaPaca_",15), 
                                                            c(rep("2D_",3), rep("3DO_", 3), 
                                                              rep("3DY_", 3), rep("PDX2D_",3), 
                                                              rep("PDX3D_", 3)), 
                                                            rep(c("1", "2", "3"),5))
colnames(temp)[1:18] <- c(str_c(rep("PANC1_",12), 
                          c(rep("2D_",4), rep("3DO_", 4), 
                            rep("3DY_", 4)),
                          rep(1:4,3)),
                          str_c(rep("PANC1_",6),
                                c(rep("PDX2D_",3), 
                                rep("PDX3D_", 3)), 
                          rep(1:3,2)))

temp <- temp[,order(colnames(temp))]

#confirm
colnames(temp)

#long format
df_long <- pivot_longer(temp, cols=1:48, names_to = "Sample", values_to = "Intensity", values_drop_na = FALSE)

#Explores NAs
summary <- df_long %>% group_by(Sample) %>% summarise("NAs" = sum(is.na(Intensity)), 
                                                      "nonNAs" = sum(!is.na(Intensity)),
                                                      "%Missing" = NAs/nonNAs*100)
#correlation matrixes
corrplot::corrplot(cor(log2(drop_na(temp))),
                   method = "color",
                   hclust.method = "ward.D",
                   is.corr = FALSE)

#Heatmap global
pheatmap::pheatmap(log2(drop_na(temp)), 
                   scale = "row", 
                   clustering_distance_cols = "correlation",
                   show_rownames = FALSE)

# boxplots of tech and biol rep
boxplot(log2(df[-1]), las=2)
boxplot(log2(temp),las=2)


df_wide <- log2(temp) %>% drop_na()
# Prepare matrix
mat <- as.matrix(df_wide)
