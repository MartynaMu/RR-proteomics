# 2D VS 3D YOUNG VS 3D OLD VS PDX - collaboration with RR
# Downstream processing of DIANN output
# INPUT: DIANN quant - FASTA with shared HUMAN and MOUSE proteins 
# Martyna Muszczek, 05.12.23

#Libraries
library(tidyverse)
library(diann)
library(smartPep)

#Import
pg_matrix <- read.delim("data/report.pg_matrix.tsv")

#store protein info in a separate df
protein_info <- pg_matrix[1:5]
df <- pg_matrix[c(3, 6:41)]

#tidy up colnames
temp <- merge2reps(df)
colnames(temp)[-1] <- str_sub(colnames(df)[-1], start=38) %>% 
  str_replace_all(pattern = "_1.+$|_2.+$|_2_repeat.+$", replacement = "") %>% 
  unique()

#long format
df_long <- pivot_longer(temp, cols=2:19, names_to = "Sample", values_to = "Intensity", values_drop_na = FALSE)

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
pheatmap::pheatmap(log2(drop_na(temp[-1])), 
                   scale = "row", 
                   clustering_distance_cols = "correlation",
                   show_rownames = FALSE)

#*Compare report pg and report.tsv after rpackage logarithm MaxLFQ


#
boxplot(log2(df[-1]))
boxplot(log2(temp[-1]))

#how many proteins are of mouse origin in PDX?
#Explores NAs
summary2 <- df_long %>% 
  group_by(Sample) %>% 
  drop_na() %>%
  summarise("NrofMurine" = sum(str_detect(Protein.Names,
                                        pattern = "_MOUSE")),
            "NrofHuman" = sum(str_detect(Protein.Names,
                                        pattern = "_HUMAN")),
            "%ofMurine" = sum(NrofMurine/(NrofHuman+NrofMurine)*100),
            "All" = NrofMurine+NrofHuman)

#filter protein groups
df_no_pg <- filter(temp, str_detect(temp$Protein.Names, ";") == FALSE)

summary2 <- pivot_longer(df_no_pg, cols=2:19, names_to = "Sample", values_to = "Intensity") %>% 
  group_by(Sample) %>% 
  drop_na() %>%
  summarise("NrofMurine" = sum(str_detect(Protein.Names,
                                          pattern = "_MOUSE")),
            "NrofHuman" = sum(str_detect(Protein.Names,
                                         pattern = "_HUMAN")),
            "%ofMurine" = sum(NrofMurine/(NrofHuman+NrofMurine)*100))

