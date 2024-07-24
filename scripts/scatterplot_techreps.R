#pairsplot for tech reps
library("tidyverse")
library("ggplot2")
library("ggpubr")

colnames(df)[-1] <- str_sub(colnames(df)[-1], start=37) %>% 
  str_replace_all(pattern = ".mzML.dia|_repeat.mzML.dia", replacement = "") %>%
  str_replace_all(pattern = "^_", replacement = "s")

theme_set(theme_minimal())

ggplot(log2(df[-1]))+
  geom_point(aes(x=s2D1_1, y=s2D1_2), color = "blue", alpha = .5)+
  geom_abline(color = "red", size = 1)

plot <- ggscatter(log2(df[-1]), x="s2D1_1", y="s2D1_2",
                add = "reg.line",
                add.params = list(color = "blue"),
                color = "darkgray")
plot + stat_cor(method = "pearson", p.accuracy = 0.001, r.accuracy = 0.01)

tech_pairsplot <- function(wide_df, cols, isLogTransformed = FALSE) {
  require("ggpubr")
  wide_df[cols] <- log2(wide_df[cols])
  #empty list to append new plots to
  scatter_list <- list()
  #iterate every 2 techreps 
  for (i in seq.int(first(cols), last(cols), 2)) {
    p1 <- ggscatter(wide_df[cols],
                        x = wide_df[[i]],
                        y = wide_df[[i+1]],
                        add = "reg.line",
                        add.params = list(color = "blue"), 
                        color = "darkgray")+
                stat_cor(method = "pearson", p.accuracy = 0.001, r.accuracy = 0.01)
    #and add the plot into list
    scatter_list[[length(scatter_list)+1]] <- p1
  }
  return(scatter_list)
}
temp <- tech_pairsplot(df, 2:37)

ggarrange(plotlist = temp)