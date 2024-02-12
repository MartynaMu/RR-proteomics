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
groups <- list(1:4,5:8,9:12,13:15,16:18)
names(groups) <- colnames(design)

# Generalized
groups2 <- list(1:4,5:12,13:18)
names(groups2) <- colnames(design2)

#take 2 experimental groups from matrix
x <- rowMeans(qnorm[,c(groups$x3D_old)])
y <- rowMeans(qnorm[,c(groups$PDX_3D)])

M <- x - y
A <- (x + y)/2

#plot and save
png("figures/maplot.png")
maplot <- ma.plot(A = A, M = M, cex = 1)
title("... normalisation")
dev.off()

#blue line - theoretical perfect loess curve
#red line - actual loess curve
