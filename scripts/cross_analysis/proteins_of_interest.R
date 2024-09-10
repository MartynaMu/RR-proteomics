# Boxplots of proteins of interest with statistics

library(tidyverse)
library(ggpubr)

mat.long <- mat %>% rownames_to_column("Gene") %>%
  pivot_longer(2:49,
    names_to = c("Line", "Sample"),
    values_to = "Intensity",
    names_sep = "(?<=PANC|MIAPACA|CFPAC)\\_")

mat.long <- mat.long %>% separate(col = "Sample", into = c("Group", "Rep"), sep = "_(?=[[:digit:]])")

boxplot.comp <- list(c("2D", "3D_YOUNG"),
                     c("3D_YOUNG", "3D_OLD"),
                     c("3D_OLD", "2D_XENO"),
                     c("2D_XENO", "3D_XENO"),
                     c("2D", "3D_OLD"), 
                     c("3D_OLD", "3D_XENO"),
                     c("2D", "2D_XENO"),
                     c("2D", "3D_XENO"))

curr.gene = "PODXL"
p <- mat.long %>% 
  filter(Gene == curr.gene) %>% 
  ggline("Group", 
     "Intensity", 
     color = "Line",
     add=c("mean_se","dotplot"),
     title = paste0(curr.gene," expression")) |>
  ggpar(legend = "none", x.text.angle = 45) |>
  facet(facet.by = "Line", panel.labs.background = list(color="white", fill="white"))+
  stat_compare_means(method = "anova", label.y.npc = "bottom")+
  stat_compare_means(method = "t.test", comparisons = boxplot.comp, label = "p.signif")+
  font("x.text", size=8)
ggsave(p, filename = paste0(curr.gene, ".png"), device = "png",path = "figures/allruns/final_quant", width = 8, height = 5, units = "in", scale = 1, dpi=100)
p  

