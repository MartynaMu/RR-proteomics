# x2D_x3Dyoung = x2D - x3D_young, #1
# x2D_x3Dold = x2D - x3D_old, #2
# x2D_PDX2D = x2D - PDX_2D, #3
# x3Dyoung_x3Dold=x3D_young - x3D_old, #4
# x3Dold_PDX3D = x3D_old - PDX_3D, #5
# x3Dyoung_xPDX3D = x3D_young - PDX_3D, #6

# Libs-------------------------------------------------------------------------
library(tidyverse)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggrepel)

coefs <- seq.int(1,6,1)
names(coefs) <- c("2Dv3Dy", "2Dv3Do", "2DvPDX2D", "3Dyv3Do", "3DovPDX3D", "3DyvPDX3D")

for (i in coefs) {
  if (dir.exists("figures/volcanoes/") == FALSE) {
    dir.create("figures/volcanoes/",showWarnings = TRUE,recursive = TRUE)
  }
  if (dir.exists("data/comparison_stats/") == FALSE) {
    dir.create("data/comparison_stats/",showWarnings = TRUE,recursive = TRUE)
  }

  deps <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(qnorm))
  
  volc_data <- deps
  volc_data$Category <- "Not significant"
  volc_data$Category[volc_data$logFC <= -1.3 & volc_data$adj.P.Val <= 0.05] <- "Down-regulated"
  volc_data$Category[volc_data$logFC >= 1.3 & volc_data$adj.P.Val <= 0.05] <- "Up-regulated"
  volc_data$Label <- NA
  volc_data$Label[volc_data$Category != "Not significant"] <- volc_data$ID[volc_data$Category != "Not significant"]
  
  upreg <- volc_data$ID[volc_data$logFC >= 1.3 & volc_data$adj.P.Val <= 0.05]
  downreg <- volc_data$ID[volc_data$logFC <= -1.3 & volc_data$adj.P.Val <= 0.05]
  
  write.table(file = paste0("data/comparison_stats/", names(coefs[i]), ".tsv"), 
              x = volc_data, 
              quote = FALSE, 
              sep = "\t",
              row.names = FALSE)
  
  theme_set(theme_light())
  v <- ggplot(volc_data,
              aes(x = logFC,
                  y = (-1*log10(adj.P.Val)),
                  color = Category,
                  label = Label))+
    geom_point(alpha = .5)+
    scale_color_manual(values = c("Down-regulated" = "blue", "Not significant" = "gray", "Up-regulated" = "red"))+
    labs(title = names(coefs[i]),
         y = "-log10 adj. p-value",
         x = "log2 fold change")+
    geom_hline(yintercept = 1.3, lty=2)+
    geom_vline(xintercept = c(-1.3,1.3), lty=2)+
    scale_x_continuous(limits = c(-5.5,5.5))+
    scale_y_continuous(limits = c(0,10.5))
  v+
    theme(legend.position = "top",
          text = element_text(size=15))+
    geom_text_repel(colour = "black", max.overlaps = 15, size = 3.5)+
    annotate("text", x = -5.5, y=0, col= "blue",label=paste0("n=",length(downreg)))+
    annotate("text", x = 5.5, y=0, col= "red",label=paste0("n=",length(upreg)))
  
  ggsave(filename = paste0("volcano", names(coefs[i]), ".png"),
         path = "figures/volcanoes/",
         device = "png",
         width = 640,
         height = 660,
         units = "px",
         dpi = 100)
}
