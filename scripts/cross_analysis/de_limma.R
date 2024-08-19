library(limma)
mat <- as.matrix(qnorm)

design <- model.matrix(~ 0+factor(c(rep(1,4),
                                    rep(2,4),
                                    rep(3,4),
                                    rep(4,3),
                                    rep(5,3),
                                    rep(6,3),
                                    rep(7,3),
                                    rep(8,3),
                                    rep(9,3),
                                    rep(10,3),
                                    rep(11,3),
                                    rep(12,3),
                                    rep(13,3),
                                    rep(14,3),
                                    rep(15,3))))

colnames(design) <- c("PANC.x2D", "PANC.x3D_young", "PANC.x3D_old", "PANC.PDX_2D", "PANC.PDX_3D",
                      "MIAPACA.x2D", "MIAPACA.x3D_young", "MIAPACA.x3D_old", "MIAPACA.PDX_2D", "MIAPACA.PDX_3D",
                      "CFPAC.x2D", "CFPAC.x3D_young", "CFPAC.x3D_old", "CFPAC.PDX_2D", "CFPAC.PDX_3D")

row.names(design) <- colnames(mat)

fit <- lmFit(mat, design)
cont.matrix <- makeContrasts(PANC.x2D_x3Dyoung = PANC.x3D_young - PANC.x2D, #1
                             PANC.x2D_x3Dold = PANC.x3D_old - PANC.x2D, #2
                             PANC.x2D_PDX2D = PANC.PDX_2D - PANC.x2D, #3
                             PANC.x3Dyoung_x3Dold = PANC.x3D_old - PANC.x3D_young, #4
                             PANC.x3Dold_PDX3D = PANC.PDX_3D - PANC.x3D_old, #5
                             PANC.x3Dyoung_xPDX3D = PANC.PDX_3D - PANC.x3D_young, #6
                             MIAPACA.x2D_x3Dyoung = MIAPACA.x3D_young - MIAPACA.x2D, #7
                             MIAPACA.x2D_x3Dold = MIAPACA.x3D_old - MIAPACA.x2D, #8
                             MIAPACA.x2D_PDX2D = MIAPACA.PDX_2D - MIAPACA.x2D, #9
                             MIAPACA.x3Dyoung_x3Dold = MIAPACA.x3D_old - MIAPACA.x3D_young, #10
                             MIAPACA.x3Dold_PDX3D = MIAPACA.PDX_3D - MIAPACA.x3D_old, #11
                             MIAPACA.x3Dyoung_xPDX3D = MIAPACA.PDX_3D - MIAPACA.x3D_young, #12
                             CFPAC.x2D_x3Dyoung = CFPAC.x3D_young - CFPAC.x2D, #13
                             CFPAC.x2D_x3Dold = CFPAC.x3D_old - CFPAC.x2D, #14
                             CFPAC.x2D_PDX2D = CFPAC.PDX_2D - CFPAC.x2D, #15
                             CFPAC.x3Dyoung_x3Dold = CFPAC.x3D_old - CFPAC.x3D_young, #16
                             CFPAC.x3Dold_PDX3D = CFPAC.PDX_3D - CFPAC.x3D_old, #17
                             CFPAC.x3Dyoung_xPDX3D = CFPAC.PDX_3D - CFPAC.x3D_young, #18
                             levels = design)

fit2 <- contrasts.fit(fit, cont.matrix)

fit_bayes <- eBayes(fit2) 

topTable(fit_bayes, coef = 1)
