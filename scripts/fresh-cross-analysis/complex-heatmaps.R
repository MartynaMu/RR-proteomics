library(ComplexHeatmap)

ha1 <- HeatmapAnnotation(df = annot_panc, show_annotation_name = FALSE, 
                         show_legend = c(Condition.L1 = FALSE, Condition.L2 = FALSE),
                         gp = gpar(col="black"),
                         col = list(Condition.L1 = c("2D" = "#FB9A99", "3D.young" = "#B2DF8A", 
                                    "3D.old" = "#33A02C", "PDX.2D" = "#A6CEE3", "PDX.3D" = "#1F78B4"),
                                    Condition.L2 = c("2D" = "#FD5855", "3D" = "#337F2E", "PDX" = "#2E6080")))


ht1 <- df.panc %>% 
  filter(rownames(df.panc) %in% term[[1]]) %>%
  t() %>% 
  scale() %>% 
  t() %>% 
  as.data.frame() %>%
  drop_na() %>% 
  Heatmap(top_annotation = ha1,
          show_row_names = TRUE,
          show_column_names = FALSE,column_order = order(colnames(df.panc)),
          show_row_dend = FALSE,
          cluster_columns = FALSE,
          column_split = annot_col_panc$Condition.L1,
          column_title = NULL,
          border = TRUE, 
          rect_gp = gpar(color="black"),
          heatmap_legend_param = gpar(title = NULL),
          row_names_gp = gpar(fontsize = 7))

svg("figures/cross-analysis/detailed-heatmaps/app-panc.svg", width = 5, height = 5)
draw(ht1)
dev.off()


ha2 <- HeatmapAnnotation(df = annot_miapaca, show_annotation_name = FALSE, 
                         show_legend = c(Condition.L1 = FALSE, Condition.L2 = FALSE),
                         gp = gpar(col="black"),
                         col = list(Condition.L1 = c("2D" = "#FB9A99", "3D.young" = "#B2DF8A", 
                                    "3D.old" = "#33A02C", "PDX.2D" = "#A6CEE3", "PDX.3D" = "#1F78B4"),
                                    Condition.L2 = c("2D" = "#FD5855", "3D" = "#337F2E", "PDX" = "#2E6080")))
ht2 <- df.miapaca %>% 
  filter(rownames(df.miapaca) %in% term[[1]]) %>%
  t() %>% 
  scale() %>% 
  t() %>% 
  as.data.frame() %>%
  drop_na() %>% 
  Heatmap(top_annotation = ha2,
          show_row_names = TRUE,
          show_column_names = FALSE,
          show_row_dend = FALSE,
          cluster_columns = FALSE,
          column_split = annot_miapaca$Condition.L1,
          column_title = NULL,
          border = TRUE, 
          rect_gp = gpar(color="black"),
          heatmap_legend_param = gpar(title = NULL))

svg("figures/cross-analysis/detailed-heatmaps/app-miapaca.svg", width = 5, height = 5)
draw(ht2)
dev.off()

ha3 <- HeatmapAnnotation(df = annot_cfpac, show_annotation_name = FALSE, 
                         show_legend = c(Condition.L1 = FALSE, Condition.L2 = FALSE),
                         gp = gpar(col="black"),
                         col = list(Condition.L1 = c("2D" = "#FB9A99", "3D.young" = "#B2DF8A", 
                                    "3D.old" = "#33A02C", "PDX.2D" = "#A6CEE3", "PDX.3D" = "#1F78B4"),
                                    Condition.L2 = c("2D" = "#FD5855", "3D" = "#337F2E", "PDX" = "#2E6080")))
ht3 <- df.cfpac %>% 
  filter(rownames(df.cfpac) %in% term[[1]]) %>%
  t() %>% 
  scale() %>% 
  t() %>% 
  as.data.frame() %>%
  drop_na() %>% 
  Heatmap(top_annotation = ha3,
          show_row_names = TRUE,
          show_column_names = FALSE, 
          show_row_dend = FALSE,
          cluster_columns = FALSE,
          column_split = annot_cfpac$Condition.L1,
          column_title = NULL,
          border = TRUE, 
          rect_gp = gpar(color="black"),
          heatmap_legend_param = gpar(title = NULL))

svg("figures/cross-analysis/detailed-heatmaps/app-cfpac.svg", width = 5, height = 6)
draw(ht3)
dev.off()


# arrange
ht1 <- grid.grabExpr(draw(ht1))
ht2 <- grid.grabExpr(draw(ht2))
ht3 <- grid.grabExpr(draw(ht3))
ggarrange(ht1,ht2,ht3, ncol = 3)
