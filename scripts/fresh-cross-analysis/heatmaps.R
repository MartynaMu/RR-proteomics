df.panc <- mutate(df.panc, Genes = str_replace_all(Genes, pattern = "_HUMAN", replacement = ""))
df.panc <- column_to_rownames(df.panc, var = "Genes")

df.miapaca <- mutate(df.miapaca, Genes = str_replace_all(Genes, pattern = "_HUMAN", replacement = ""))
df.miapaca <- column_to_rownames(df.miapaca, var = "Genes")

df.cfpac <- mutate(df.cfpac, Genes = str_replace_all(Genes, pattern = "_HUMAN", replacement = ""))
df.cfpac <- column_to_rownames(df.cfpac, var = "Genes")

p <- df.panc %>%
      drop_na() %>%
      pheatmap::pheatmap(scale = "row",
                         annotation_col = annot_panc,
                         annotation_colors = annot_colors,
                         annotation_names_col = FALSE,
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         col=colorRampPalette(c("blue","white","red"))(100),
                         border_color = "black",
                         main = "PANC-1 global hm")
ggsave(plot = p$gtable, "figures/cross-analysis/panc-hm.svg", device = "svg", scale = 1, dpi = 100, width=7, height=7)

p <- df.miapaca %>%
  drop_na() %>%
  pheatmap::pheatmap(scale = "row",
                     annotation_col = annot_miapaca,
                     annotation_colors = annot_colors,
                     annotation_names_col = FALSE,
                     show_rownames = FALSE,
                     show_colnames = FALSE,
                     col=colorRampPalette(c("blue","white","red"))(100),
                     border_color = "black",
                     main = "MIA-PaCa-2 global hm")
ggsave(plot = p$gtable, "figures/cross-analysis/miapaca-hm.svg", device = "svg", scale = 1, dpi = 100, width=7, height=7)

p <- df.cfpac %>%
  drop_na() %>%
  pheatmap::pheatmap(scale = "row",
                     annotation_col = annot_cfpac,
                     annotation_colors = annot_colors,
                     annotation_names_col = FALSE,
                     show_rownames = FALSE,
                     show_colnames = FALSE,
                     col=colorRampPalette(c("blue","white","red"))(100),
                     border_color = "black",
                     main = "CFPAC-1 global hm")
ggsave(plot = p$gtable, "figures/cross-analysis/cfpac-hm.svg", device = "svg", scale = 1, dpi = 100, width=7, height=7)


# gene sets heatmaps

term <- filter(term.desc.bp, ID == "0019882") %>% pull(Genes)
names(term) <- "0019882"

p <- df.panc %>%
  filter(rownames(df.panc) %in% term[[1]]) %>%
  pheatmap::pheatmap(scale = "row",
                     annotation_col = annot_panc,
                     annotation_colors = annot_colors,
                     annotation_names_col = FALSE,
                     show_rownames = TRUE,
                     show_colnames = FALSE,
                     col=colorRampPalette(c("blue","white","red"))(100),
                     border_color = "black",
                     main = "PANC-1 global hm")
ggsave(plot = p$gtable, "figures/cross-analysis/panc-hm.svg", device = "svg", scale = 1, dpi = 100, width=7, height=7)

