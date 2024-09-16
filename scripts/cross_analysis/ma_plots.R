library(tidyverse)
library(affy)
library(ggpubr)

panc.gr <- list(1:4,5:8,9:12,13:15,16:18)
names(panc.gr) <- annot_col[1:18,1] %>% unique()
miapaca.gr <- list(19:21, 22:24, 25:27, 28:30, 31:33)
names(miapaca.gr) <- annot_col[19:33,1] %>% unique()
cfpac.gr <- list(34:36, 37:39, 40:42, 43:45, 46:48)
names(cfpac.gr) <- annot_col[34:48,1] %>% unique()

ds = zscored
is.norm = "Medians-centered-cell-lines-zscored"
comb <- expand.grid(1:5,1:5) |> filter(Var1 != Var2) # every possible combination of groups comparison
titles = c("PANC1", "MIAPACA", "CFPAC")
caption = "M - mean expression in two groups; A - log2FC between groups; blue - perfect loess curve; red - actual loess curve"

ma_plots <- list()
counter = 0
for(l in list(panc.gr, miapaca.gr, cfpac.gr)) {
  counter <<- counter + 1
  for (i in 1:nrow(comb)) {
    ind1 <- l[[comb[i,1]]] # create indices of current groups in each cell line
    ind2 <- l[[comb[i,2]]]
    x <- rowMeans(ds[,ind1])
    y <- rowMeans(ds[,ind2])
    M <- x - y
    A <- (x + y)/2
    temp <- data.frame("M" = M, "A" = A)
    title <- paste0(names(l[comb[i,1]])," - ", names(l[comb[i,2]]))
    maplot <- ggscatter(temp, x = "A", y = "M",
              shape = 1,
              size = 1,
              add = "loess", 
              add.params = list(color = "red", size=.5),
              title = title)+
      geom_hline(yintercept = 0, colour = "blue", lty=1, linewidth = .5)+
      font("axis.title", size=8)+
      font("xy.text", size=8)+
      font("title", size = 10)
    ma_plots <- c(ma_plots, list(maplot))
  }
  plot <- do.call(ggarrange, c(ma_plots, ncol=4, nrow=5)) |>
    annotate_figure(top = text_grob(paste(titles[counter], is.norm, sep = " - "), size = 15),
                    bottom = text_grob(caption, size = 8, just="left", x=0.01))
  filename = paste(titles[counter],is.norm,"maplot.png", sep="-")
  ggsave(plot, filename = filename, device = "png", path = "figures/allruns/final_quant/normalization/", bg = "white", scale = 1, width = 10, height = 8, dpi = 100)
  ma_plots <- list() # empty list to store plots in
}

#blue line - theoretical perfect loess curve
#red line - actual loess curve


