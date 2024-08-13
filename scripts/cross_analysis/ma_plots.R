panc.gr <- list(1:4,5:8,9:12,13:15,16:18)
names(panc.gr) <- annot_col[1:18,1] %>% unique()
miapaca.gr <- list(19:21, 22:24, 25:27, 28:30, 31:33)
names(miapaca.gr) <- annot_col[19:33,1] %>% unique()
cfpac.gr <- list(34:36, 37:39, 40:42, 43:45, 46:48)
names(cfpac.gr) <- annot_col[34:48,1] %>% unique()


#take 2 experimental groups from matrix
x <- rowMeans(test[,c(cfpac.gr$`2D`)])
y <- rowMeans(test[,c(panc.gr$`2D`)])

x <- rowMeans(qnorm[,c(cfpac.gr$`2D`)])
y <- rowMeans(qnorm[,c(panc.gr$`2D`)])

M <- x - y
A <- (x + y)/2

#plot and save
library(affy)
maplot <- ma.plot(A = A, M = M, cex = 1)

#blue line - theoretical perfect loess curve
#red line - actual loess curve
