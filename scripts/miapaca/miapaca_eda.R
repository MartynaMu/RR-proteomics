
# Norm ---------------------------------------------------------------------
qnorm <- limma::normalizeQuantiles(drop_na(int_mat))
qnorm <- limma::normalizeCyclicLoess(drop_na(int_mat))
qnorm <- limma::normalizeMedianValues(drop_na(int_mat))

# Correlation----------------------------------------------------------
## Pre-norm=============
int_mat %>% 
  drop_na() %>% 
  log2() %>% 
  cor() %>%
  corrplot::corrplot(method = "color",
                   hclust.method = "ward.D",
                   is.corr = FALSE)

## Post-norm============
qnorm %>% 
  log2() %>% 
  cor() %>%
  corrplot::corrplot(method = "color",
                     hclust.method = "ward.D",
                     is.corr = FALSE)


# PCA -------------------------------------------------------------------
library(PCAtools)
## Pre-norm ===================================================
p <- pca(drop_na(int_mat), scale = TRUE)

biplot(p, 
       lab = colnames(int_mat), 
       # labSize = 5, 
      # shape = "Condition.L3",
      # shapekey = c("CFPAC" = 15, "MiaPaca" = 17, "PANC1" = 8),
       showLoadings = FALSE,
       boxedLoadingsNames = FALSE, 
      # colby = "Condition.L1",
      # colkey = annot_colors$Condition.L1,
       title = "Pre-normalisation",
       # encircle = TRUE,
       legendPosition = "right",
       hline = 0, 
       vline = 0)

## Post-norm ======================================
p <- pca(qnorm, scale = TRUE)
biplot(p, 
       lab = colnames(qnorm), 
       # labSize = 5, 
       # shape = "Condition.L3",
       # shapekey = c("CFPAC" = 15, "MiaPaca" = 17, "PANC1" = 8),
       showLoadings = FALSE,
       boxedLoadingsNames = FALSE, 
       # colby = "Condition.L1",
       # colkey = annot_colors$Condition.L1,
       title = "Post-normalisation",
       # encircle = TRUE,
       legendPosition = "right",
       hline = 0, 
       vline = 0)

# MA plot ---------------------
# M = log2(x/y) aka log2FC
# A = (log2(x) + log2(y))/2 = log2(xy)*1/2, aka means of means=average expression
# where x and y are respectively the mean of the two groups being compared.
library(affy)

# x2D_x3Dyoung = x2D - x3D_young, #1
# x2D_x3Dold = x2D - x3D_old, #2
# x2D_PDX2D = x2D - PDX_2D, #3
# x3Dyoung_x3Dold=x3D_young - x3D_old, #4
# x3Dold_PDX3D = x3D_old - PDX_3D, #5
# x3Dyoung_xPDX3D = x3D_young - PDX_3D, #6

# In depth
groups <- list(1:3,4:6,7:9,10:12, 13:15)
names(groups) <- colnames(design)

# Generalized
# groups2 <- list(1:4,5:12,13:18)
# names(groups2) <- colnames(design2)

## Pre-norm ===============
#take 2 experimental groups from matrix
x <- rowMeans(drop_na(int_mat[,c(groups$PDX_2D)]))
y <- rowMeans(drop_na(int_mat[,c(groups$PDX_3D)]))

M <- x - y
A <- (x + y)/2

#plot and save
png("figures/miapaca/maplot.png")
maplot <- ma.plot(A = A, M = M, cex = 1)
title("Pre normalisation")
dev.off()

#blue line - theoretical perfect loess curve
#red line - actual loess curve

## Post-norm =====
#take 2 experimental groups from matrix
x <- rowMeans(qnorm[,c(groups$PDX_2D)])
y <- rowMeans(qnorm[,c(groups$PDX_3D)])

M <- x - y
A <- (x + y)/2

#plot and save
png("figures/miapaca/maplot2.png")
maplot <- ma.plot(A = A, M = M, cex = 1)
title("Post normalisation")
dev.off()
