# Libs-------------------------------------------------------------------------
library(tidyverse)
library(ggrepel)
# MSqrob/Hurdle fc--------------------------------------------------------------
coefs <- seq.int(1,ncol(fc_msqrob),2)
names(coefs) <- fc_msqrob %>% dplyr::select(matches("FC")) %>% colnames() %>% str_replace_all("_FC", "")

## set plot axis limits-----------------------------------------------------------
max.y <- fc_msqrob %>% dplyr::select(ends_with("qval")) %>% drop_na() %>% summarise_all(.funs = min) %>% min() %>% log10()*(-1)+0.5
min.x <- fc_msqrob %>% dplyr::select(ends_with("FC")) %>% drop_na() %>% summarise_all(.funs = min) %>% min()-0.5
max.x <- fc_msqrob %>% dplyr::select(ends_with("FC")) %>% drop_na() %>% summarise_all(.funs = max) %>% max()+0.5

## loop---------------------------------------------------------------------------
for (i in coefs) {
  if (dir.exists("figures/panc1/hurdle_msqrob_fc/volcanoes/") == FALSE) {
    dir.create("figures/panc1/hurdle_msqrob_fc/volcanoes/",showWarnings = TRUE,recursive = TRUE)
  }
  if (dir.exists("data/panc1/hurdle_msqrob_fc/comparison_stats/") == FALSE) {
    dir.create("data/panc1/hurdle_msqrob_fc/comparison_stats/",showWarnings = TRUE,recursive = TRUE)
  }
  
  message(sprintf("Preparing %s volcano plot", names(coefs[coefs==i])))
  
  volc_data <- fc_msqrob[c(i, i+1)]
  colnames(volc_data) <- c("logFC", "adj.P.Val")
  volc_data <- rownames_to_column(volc_data, var = "ID")
  volc_data$Category <- "Not significant"
  volc_data$Category[volc_data$logFC <= -1.3 & volc_data$adj.P.Val <= 0.05] <- "Down-regulated"
  volc_data$Category[volc_data$logFC >= 1.3 & volc_data$adj.P.Val <= 0.05] <- "Up-regulated"
  volc_data$Label <- NA
  volc_data$Label[volc_data$Category != "Not significant"] <- volc_data$ID[volc_data$Category != "Not significant"]
  
  upreg <- volc_data$ID[volc_data$Category == "Up-regulated"]
  downreg <- volc_data$ID[volc_data$Category == "Down-regulated"]
  
  write.table(file = paste0("data/panc1/hurdle_msqrob_fc/comparison_stats/", names(coefs[coefs==i]), ".tsv"), 
              x = volc_data, 
              quote = FALSE, 
              sep = "\t",
              row.names = FALSE)
  
  #set comparison names
  comp <- str_split_1(names(coefs[coefs==i]), pattern = "_vs_")

  theme_set(theme_light())
  v <- ggplot(volc_data,
              aes(x = logFC,
                  y = (-1*log10(adj.P.Val)),
                  color = Category,
                  label = Label))+
    geom_point(alpha = .5)+
    scale_color_manual(values = c("Down-regulated" = "blue", "Not significant" = "gray", "Up-regulated" = "red"))+
    labs(title = names(coefs[coefs==i]),
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
  
  ggsave(filename = paste0("volcano_", names(coefs[coefs==i]), ".png"),
         path = "figures/panc1/hurdle_msqrob_fc/volcanoes/",
         device = "png",
         width = 640,
         height = 660,
         units = "px",
         dpi = 100)
}


# Int matrix fc (limma) ----------------------------------------------------------------
coefs <- seq.int(1,6,1)
names(coefs) <- c("2Dv3Dy", "2Dv3Do", "2DvPDX2D", "3Dyv3Do", "3DovPDX3D", "3DyvPDX3D")

## set plot axis limits---------------------------------------------------------
# min.x <- c()
# for (i in coefs) {
#   x <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(int_mat)) %>% 
#     dplyr::select("logFC") %>% drop_na() %>%
#     min()
#   min.x <- c(min.x, x)
# }
# min.x <- min(min.x)-0.5
# 
# max.x <- c()
# for (i in coefs) {
#   x <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(int_mat)) %>% 
#     dplyr::select("logFC") %>% drop_na() %>%
#     max()
#   max.x <- c(max.x, x)
# }
# max.x <- max(max.x)+0.5
# 
# max.y <- c()
# for (i in coefs) {
#   y <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(int_mat)) %>% 
#     dplyr::select("adj.P.Val") %>% drop_na() %>%
#     min()
#   max.y <- c(max.y, y)
# }
# max.y <- min(max.y) %>% log10()*(-1)+0.5

max.x <- 7.5
min.x <- -7.5
max.y <- 10.5

## loop -----------------------------------------------------------------------
for (i in coefs) {
  if (dir.exists("figures/panc1/int_mat_fc/volcanoes/") == FALSE) {
    dir.create("figures/panc1/int_mat_fc/volcanoes/",showWarnings = TRUE,recursive = TRUE)
  }
  if (dir.exists("data/panc1/int_mat_fc/comparison_stats/") == FALSE) {
    dir.create("data/panc1/int_mat_fc/comparison_stats/",showWarnings = TRUE,recursive = TRUE)
  }
  
  deps <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(int_mat))
  
  volc_data <- deps
  volc_data$Category <- "Not significant"
  volc_data$Category[volc_data$logFC <= -1.3 & volc_data$adj.P.Val <= 0.05] <- "Down-regulated"
  volc_data$Category[volc_data$logFC >= 1.3 & volc_data$adj.P.Val <= 0.05] <- "Up-regulated"
  volc_data$Label <- NA
  volc_data$Label[volc_data$Category != "Not significant"] <- volc_data$ID[volc_data$Category != "Not significant"]
  
  upreg <- volc_data$ID[volc_data$Category == "Up-regulated"]
  downreg <- volc_data$ID[volc_data$Category == "Down-regulated"]
  
  write.table(file = paste0("data/panc1/int_mat_fc/comparison_stats/", names(coefs[i]), ".tsv"), 
              x = volc_data, 
              quote = FALSE, 
              sep = "\t",
              row.names = FALSE)
  
  #set comparison name
  comp <- str_split_1(names(coefs[i]), pattern = "v")
  
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
         path = "figures/panc1/int_mat_fc/volcanoes/",
         device = "png",
         width = 640,
         height = 660,
         units = "px",
         dpi = 100)
}

