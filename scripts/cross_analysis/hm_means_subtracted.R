library(tidyverse)
df <- read.delim("data/all_lines/a_Hurdle_Msqrob.tsv", header = TRUE, sep = "\t")
df <- df[1:49]

df$Genes <- str_replace_all(df$Genes, "_HUMAN", replacement = "")
df <- column_to_rownames(df, var = "Genes")
df <- limma::normalizeQuantiles(df)
df <- rownames_to_column(df, var = "Genes")

panc.ind <- colnames(df) %>% str_which(pattern="PANC")
miapaca.ind <- colnames(df) %>% str_which(pattern="MIAPACA")
cfpac.ind <- colnames(df) %>% str_which(pattern="CFPAC")

test <- mutate(rowwise(df), PANC.mean = mean(c_across(all_of(panc.ind)),na.rm=TRUE),
               MIAPACA.mean = mean(c_across(all_of(miapaca.ind)), na.rm = TRUE),
               CFPAC.mean = mean(c_across(all_of(cfpac.ind)), na.rm = TRUE))

test <- mutate(test, across(.cols = all_of(panc.ind), .fns = function(x) (x - PANC.mean)),
               across(.cols = all_of(miapaca.ind), .fns = function(x) (x - MIAPACA.mean)),
               across(.cols = all_of(cfpac.ind), .fns = function(x) (x - CFPAC.mean)))

# Complex hm with annotations
cfpac.l1 <- rep(c("2D", "3D.young", "3D.old", "PDX.2D", "PDX.3D"), each = 3)
miapaca.l1 <- rep(c("2D", "3D.young", "3D.old", "PDX.2D", "PDX.3D"), each = 3)
panc1.l1 <- c(rep(c("2D", "3D.young", "3D.old"), each = 4), rep(c("PDX.2D", "PDX.3D"), each = 3))

cfpac.l2 <- c(rep("2D", 3), rep("3D", 6), rep("PDX", 6))
miapaca.l2 <- c(rep("2D", 3), rep("3D", 6), rep("PDX", 6))
panc1.l2 <- c(rep("2D", 4), rep("3D", 8), rep("PDX", 6))

celllines <- c(rep("PANC1", 18), rep("MiaPaca", 15), rep("CFPAC", 15))


annot_col <- data.frame(
  Condition.L1 = as.factor(
    c(panc1.l1, miapaca.l1, cfpac.l1)),
  Condition.L2 = as.factor(
    c(panc1.l2, miapaca.l2, cfpac.l2)),
  Condition.L3 = as.factor(celllines))

row.names(annot_col) <- colnames(test[2:49])

library(RColorBrewer)
annot_colors <- list(Condition.L1 = rev(brewer.pal(5, "Paired")),
                     Condition.L2 = c("firebrick1", "darkgreen", "blue"),
                     Condition.L3 = c("black", "darkgrey", "lightgray"))

names(annot_colors$Condition.L1) <- levels(annot_col$Condition.L1)
names(annot_colors$Condition.L2) <- levels(annot_col$Condition.L2)
names(annot_colors$Condition.L3) <- levels(annot_col$Condition.L3)

test[test==0] <- NA
test <- test[2:49] %>% 
  drop_na()

pheatmap::pheatmap(test, scale = "row",
                   annotation_col = annot_col,
                   annotation_colors = annot_colors,
                   main = "Global heatmap qnorm, cell line means subtracted")

