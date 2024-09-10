# SCATTERPLOT OF FC OF DEPS BETWEEN CELL LINES IN EACH COMPARISON
scatter <- list()
for (i in curr.comp) {
  temp <- FC %>% select(contains(i))
  a = colnames(temp[1])
  b = colnames(temp[2])
  c = colnames(temp[3])
  
  p <- temp %>% ggscatter(x=a, 
                          y=b,
                          shape=1,
                          size=1,
                     cor.coef = TRUE, 
                     cor.coeff.args = list(method = "pearson"))+
    font("axis.title", size=8)+
    font("xy.text", size=8)
  scatter <- c(scatter, list(p))
  p <- temp %>% ggscatter(x=a, 
                          y=c,
                          shape=1,
                          size=1,
                          cor.coef = TRUE, 
                          cor.coeff.args = list(method = "pearson"))+
    font("axis.title", size=8)+
    font("xy.text", size=8)
  scatter <- c(scatter, list(p))
  p <- temp %>% ggscatter(x=b, 
                          y=c,
                          shape=1,
                          size=1,font.label = 6,
                          cor.coef = TRUE, 
                          cor.coeff.args = list(method = "pearson"))+
    font("axis.title", size=8)+
    font("xy.text", size=8)
  scatter <- c(scatter, list(p))
}

multi.scatter <- do.call(ggarrange,c(scatter, ncol=3, nrow=6))

ggsave(multi.scatter, filename="figures/allruns/final_quant/multi.scatter.png", device = "png", dpi=150, width=8, height=10, scale=1)
