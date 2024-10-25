library(tidyverse)

# read in all quantity matrices obtained from separate DIANN runs for each line
df.panc <- read.delim("data/panc1/a_Hurdle_Msqrob.tsv", header = TRUE, sep = "\t")
df.panc <- df.panc[1:19]

df.miapaca <- read.delim("data/miapaca/a_Hurdle_Msqrob.tsv", header = TRUE, sep = "\t")
df.miapaca <- df.miapaca[1:16]

df.cfpac <- read.delim("data/cfpac/a_Hurdle_Msqrob.tsv", header = TRUE, sep = "\t")
df.cfpac <- df.cfpac[1:16]

# center medians where needed
boxplot(df.panc[1:18])

df.miapaca[df.miapaca==0] <- NA
boxplot(df.miapaca[1:15])
boxplot(df.panc[1:15])

# full join by Genes
df <- full_join(df.panc, df.miapaca, by = "Genes")
df <- full_join(df, df.cfpac, by = "Genes")

# Tidy gene names
df$Genes <- str_replace_all(df$Genes, "_HUMAN", replacement = "")
df <- column_to_rownames(df, var = "Genes")
df[df==0] <- NA

boxplot(df)
# Don't drop NAs from joint matrix as it will disturb normalization in cell lines

# order groups 
order <- colnames(df.panc) %>% str_replace_all(pattern = "MIAPACA_|PANC_|CFPAC_", replacement = "") %>% unique()

temp <- c()
for (i in order) {
  temp <- c(temp, str_which(colnames(df.panc), pattern = i))
}
df.panc <- relocate(df.panc, all_of(temp))


# annotations for heatmaps
cfpac.l1 <- rep(c("2D", "3D.young", "3D.old", "PDX.2D", "PDX.3D"), each = 3)
miapaca.l1 <- rep(c("2D", "3D.young", "3D.old", "PDX.2D", "PDX.3D"), each = 3)
panc1.l1 <- c(rep(c("2D", "3D.young", "3D.old"), each = 4), rep(c("PDX.2D", "PDX.3D"), each = 3))

cfpac.l2 <- c(rep("2D", 3), rep("3D", 6), rep("PDX", 6))
miapaca.l2 <- c(rep("2D", 3), rep("3D", 6), rep("PDX", 6))
panc1.l2 <- c(rep("2D", 4), rep("3D", 8), rep("PDX", 6))

annot_panc <- data.frame(
  Condition.L1 = as.factor(panc1.l1),
  Condition.L2 = as.factor(panc1.l2))

annot_cfpac <- data.frame(
  Condition.L1 = as.factor(cfpac.l1),
  Condition.L2 = as.factor(cfpac.l2))

annot_miapaca <- data.frame(
  Condition.L1 = as.factor(miapaca.l1),
  Condition.L2 = as.factor(miapaca.l2))

rownames(annot_panc) <- colnames(df.panc[-19])
rownames(annot_cfpac) <- colnames(df.cfpac[-16])
rownames(annot_miapaca) <- colnames(df.miapaca[-16])

library(RColorBrewer)
annot_colors_panc <- list(Condition.L1 = rev(brewer.pal(5, "Paired")),
                     Condition.L2 = c("firebrick1", "darkgreen", "blue"))

names(annot_colors$Condition.L1) <- levels(annot_panc$Condition.L1)
names(annot_colors$Condition.L2) <- levels(annot_panc$Condition.L2)

# heatmaps
p1 <- df[panc.ind]  %>%
  filter(rownames(df) %in% panc.canc$gene) %>%
  drop_na() %>%
  pheatmap::pheatmap(scale = "row",
                     annotation_col = annot_panc,
                     annotation_colors = annot_colors,
                     annotation_names_col = FALSE,
                     show_rownames = TRUE,
                     show_colnames = FALSE,
                     treeheight_row = 0,
                     cluster_cols = FALSE,
                     gaps_col = c(4,8,12),
                     col=colorRampPalette(c("blue","white","red"))(100),
                     border_color = "black")

ggsave(plot = p, "figures/cross-analysis/miapaca-panc-canc-genes.svg", device = "svg", scale = 1, dpi = 100, width=5.5, height=5.5)
