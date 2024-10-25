library(PCAtools)


annot_col_panc <- filter(annot_col, Condition.L3 == "PANC1")
p <- pca(drop_na(df.panc), metadata = annot_col_panc)
biplot(p, 
       lab = colnames(df.panc),
       colkey = annot_colors$Condition.L1,
       colby = "Condition.L1",
       encircle = TRUE,
       hline = 0,
       vline = 0)

ggsave("figures/cross-analysis/pca-panc.svg", device = "svg", scale = 1, dpi = 100, width = 6, height = 6)


annot_col_miapaca <- filter(annot_col, Condition.L3 == "MiaPaca")
p <- pca(drop_na(df.miapaca), metadata = annot_col_miapaca)
biplot(p, 
       lab = colnames(df.miapaca),
       colkey = annot_colors$Condition.L1,
       colby = "Condition.L1",
       encircle = TRUE,
       hline = 0,
       vline = 0)

ggsave("figures/cross-analysis/pca-miapaca.svg", device = "svg", scale = 1, dpi = 100, width = 6, height = 6)


annot_col_cfpac <- filter(annot_col, Condition.L3 == "CFPAC")
p <- pca(drop_na(df.cfpac), metadata = annot_col_cfpac)
biplot(p, 
       lab = colnames(df.cfpac),
       colkey = annot_colors$Condition.L1,
       colby = "Condition.L1",
       encircle = TRUE,
       hline = 0,
       vline = 0)

ggsave("figures/cross-analysis/pca-cfpac.svg", device = "svg", scale = 1, dpi = 100, width = 6, height = 6)
