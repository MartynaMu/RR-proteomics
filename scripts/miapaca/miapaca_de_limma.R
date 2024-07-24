# limma instead of hurdle/msqrob ttest
library(tidyverse)
library(limma)

# In-depth comparisons-----------------------------------------------------------
##Design table creation-------------------------------------------------------------
#where first digit in rep represents the column it will be filled with, second is nr of replicates
mat <- as.matrix(qnorm)
design <- model.matrix(~ 0+factor(c(rep(1,3),
                                    rep(2,3),
                                    rep(3,3),
                                    rep(4,3),
                                    rep(5,3))))
colnames(design) <- c("x2D", "x3D_young", "x3D_old", "PDX_2D", "PDX_3D")
row.names(design) <- colnames(mat)

##Limma fit and matrix contrasts-------------------------------------------------
fit <- lmFit(mat, design)
cont.matrix <- makeContrasts(x2D_x3Dyoung = x3D_young - x2D, #1
                             x2D_x3Dold = x3D_old - x2D, #2
                             x2D_PDX2D = PDX_2D - x2D, #3
                             x3Dyoung_x3Dold = x3D_old - x3D_young, #4
                             x3Dold_PDX3D = PDX_3D - x3D_old, #5
                             x3Dyoung_xPDX3D = PDX_3D - x3D_young, #6
                             levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)

fit_bayes <- eBayes(fit2) 


