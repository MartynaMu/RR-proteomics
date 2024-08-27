panc.gr <- list(1:4,5:8,9:12,13:15,16:18)
names(panc.gr) <- annot_col[1:18,1] %>% unique()
miapaca.gr <- list(19:21, 22:24, 25:27, 28:30, 31:33)
names(miapaca.gr) <- annot_col[19:33,1] %>% unique()
cfpac.gr <- list(34:36, 37:39, 40:42, 43:45, 46:48)
names(cfpac.gr) <- annot_col[34:48,1] %>% unique()


#take 2 experimental groups from matrix
x <- rowMeans(df[,c(miapaca.gr$`3D.young`)]) # nothing
y <- rowMeans(df[,c(miapaca.gr$`PDX.2D`)])

x <- rowMeans(med.norm[,c(miapaca.gr$`3D.young`)]) # medians norm
y <- rowMeans(med.norm[,c(miapaca.gr$`PDX.2D`)])

x <- rowMeans(test2[,c(miapaca.gr$`3D.young`)]) # medians norm then means subtracted
y <- rowMeans(test2[,c(panc.gr$`PDX.2D`)])

x <- rowMeans(test[,c(miapaca.gr$`3D.young`)]) #qnorm, then means subtracted
y <- rowMeans(test[,c(panc.gr$`PDX.2D`)])

x <- rowMeans(qnorm[,c(miapaca.gr$`3D.young`)]) # quantile normalized
y <- rowMeans(qnorm[,c(miapaca.gr$`PDX.2D`)])

x <- rowMeans(temp[,c(miapaca.gr$`3D.young`)]) # means subtracted, normalized abs medians
y <- rowMeans(temp[,c(miapaca.gr$`PDX.2D`)])

x <- rowMeans(temp2[,c(miapaca.gr$`3D.young`)]) # means subtracted, normalized cyclic loess
y <- rowMeans(temp2[,c(miapaca.gr$`PDX.2D`)])

x <- rowMeans(temp3[,c(miapaca.gr$`3D.young`)]) # means subtracted, normalized vsn
y <- rowMeans(temp3[,c(miapaca.gr$`PDX.2D`)])

M <- x - y
A <- (x + y)/2

#plot and save
library(affy)
maplot <- ma.plot(A = A, M = M, cex = 1)

#blue line - theoretical perfect loess curve
#red line - actual loess curve
