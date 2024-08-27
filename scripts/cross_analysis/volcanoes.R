# Libs-------------------------------------------------------------------------
library(tidyverse)
library(ggrepel)

# Int matrix fc (limma) ----------------------------------------------------------------
coefs <- seq.int(1,18,1)
names(coefs) <- colnames(fit_bayes$contrasts)

# ## set plot axis limits---------------------------------------------------------
# min.x <- c()
# for (i in coefs) {
#   x <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(mat)) %>%
#     dplyr::select("logFC") %>% drop_na() %>%
#     min()
#   cat(i,x)
#   min.x <- c(min.x, x)
# }
# min.x <- min(min.x)-0.5
# 
# max.x <- c()
# for (i in coefs) {
#   x <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(mat)) %>%
#     dplyr::select("logFC") %>% drop_na() %>%
#     max()
#   cat(i,x)
#   max.x <- c(max.x, x)
# }
# max.x <- max(max.x)+0.5
# 
# max.y <- c()
# for (i in coefs) {
#   y <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(mat)) %>%
#     dplyr::select("adj.P.Val") %>% drop_na() %>%
#     min()
#   cat(i,y)
#   max.y <- c(max.y, y)
# }
# max.y <- min(max.y) %>% log10()*(-1)+0.5

max.x <- 7
min.x <- -5.5
max.y <- 19.5

## loop -----------------------------------------------------------------------
for (i in coefs) {
  if (dir.exists("figures/allruns/final_quant/volcanoes/") == FALSE) {
    dir.create("figures/allruns/final_quant/volcanoes/",showWarnings = TRUE,recursive = TRUE)
  }
  if (dir.exists("data/all_lines/comparison_stats/") == FALSE) {
    dir.create("data/all_lines/comparison_stats/",showWarnings = TRUE,recursive = TRUE)
  }
  
  deps <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(qnorm))
  
  volc_data <- deps
  volc_data$Category <- "Not significant"
  volc_data$Category[volc_data$logFC <= -1.3 & volc_data$adj.P.Val <= 0.05] <- "Down-regulated"
  volc_data$Category[volc_data$logFC >= 1.3 & volc_data$adj.P.Val <= 0.05] <- "Up-regulated"
  volc_data$Label <- NA
  volc_data$Label[volc_data$Category != "Not significant"] <- volc_data$ID[volc_data$Category != "Not significant"]
  
  upreg <- volc_data$ID[volc_data$Category == "Up-regulated"]
  downreg <- volc_data$ID[volc_data$Category == "Down-regulated"]
  
  write.table(file = paste0("data/all_lines/comparison_stats/", names(coefs[i]), ".tsv"), 
              x = volc_data, 
              quote = FALSE, 
              sep = "\t",
              row.names = FALSE)
  
  #set comparison name
  comp <- str_split_1(names(coefs[i]), pattern = "_")
  
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
         x = sprintf("log2 fold change (%s/%s)", comp[2], comp[1]))+
    geom_hline(yintercept = 1.3, lty=2)+
    geom_vline(xintercept = c(-1.3,1.3), lty=2)+
    scale_x_continuous(limits = c(min.x,max.x))+
    scale_y_continuous(limits = c(0,max.y))
  v+
    theme(legend.position = "top",
          text = element_text(size=15))+
    geom_text_repel(colour = "black", max.overlaps = 15, size = 3.5)+
    annotate("text", x = -5.5, y=0, col= "blue",label=paste0("n=",length(downreg)))+
    annotate("text", x = 5.5, y=0, col= "red",label=paste0("n=",length(upreg)))
  
  ggsave(filename = paste0("volcano", names(coefs[i]), ".png"),
         path = "figures/allruns/final_quant/volcanoes/",
         device = "png",
         width = 640,
         height = 660,
         units = "px",
         dpi = 100)
}

