library(tidyverse)
library(limma)

# In-depth comparisons-----------------------------------------------------------
##Design table creation-------------------------------------------------------------
#where first digit in rep represents the column it will be filled with, second is nr of replicates
mat <- as.matrix(qnorm[1:18])
design <- model.matrix(~ 0+factor(c(rep(1,4),
                                    rep(2,4),
                                    rep(3,4),
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

##DE genes--------------------------------------------------------------------
#Coef = sth shows one chosen comparison, if all comparisons need to be done, eg. 3, anova will be done
#here 3d old vs young
deps <- topTable(fit_bayes, coef = 1, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(qnorm))
# filter on p.adj value
#deps <- filter(deps, adj.P.Val <= 0.05)
# find DEPs sign between comparisons
#anova <- topTable(fit_bayes, adjust = "BH", genelist = rownames(mat), resort.by = "adj.P.val", number = nrow(qnorm))
#filter out genes that are not sign expressed in either comparison
#anova_filtr <- anova$ProbeID[anova$adj.P.Val <= 0.05]
#deps <- filter(deps, ID %in% anova_filtr)

#save < 0,05 adj pvalue
clipr::write_clip(deps[deps$adj.P.Val < 0.05,c(1,2)])

#see how many significant DEPs overlap between comparisons - eg. downregulated
result <- decideTests(fit_bayes)
vennDiagram(result,include = "down",show.include = TRUE)

##Volcanos-------------------------------------------------------------------------
#quick volcano plot
volcanoplot(fit_bayes, coef=2, highlight = 10, names = rownames(mat))
abline(h = 1.2, v = c(-1.5,1.5), lty = 2)

# General comparisons-----------------------------------------------------------
##Design table creation-------------------------------------------------------------
#where first digit in rep represents the column it will be filled with, second is nr of replicates
design2 <- model.matrix(~ 0+factor(c(rep(1,4),
                                    rep(2,8),
                                    rep(3,6))))
colnames(design2) <- c("x2D", "x3D", "PDX")
row.names(design2) <- colnames(mat)

##Limma fit and matrix contrasts-------------------------------------------------
fit_g <- lmFit(mat, design2)
cont.matrix_g <- makeContrasts(x2D_x3D = x2D - x3D, #1
                             x2D_PDX = x2D - PDX, #2
                             x3D_PDX = x3D - PDX, #3
                             levels=design2)
fit_g <- contrasts.fit(fit_g, cont.matrix_g)

fit_bayes_g <- eBayes(fit_g) 

##DE genes--------------------------------------------------------------------
#Coef = sth shows one chosen comparison, if all comparisons need to be done, eg. 3, anova will be done
#here 3d old vs young
deps <- topTable(fit_bayes_g, coef = 3, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(qnorm))

#find DEPs sign between comparisons
anova <- topTable(fit_bayes_g, adjust = "BH", genelist = rownames(mat), resort.by = "adj.P.val", number = nrow(qnorm))

#save < 0,05 adj pvalue and paste into string-db class/ranks GO enrichment
# OR use clusterProfiler for GSEA
clipr::write_clip(deps[deps$adj.P.Val < 0.05,c(1,2)])

#see how many significant DEPs overlap between comparisons - eg. downregulated
result <- decideTests(fit_bayes_g)
vennDiagram(result,include = "down",show.include = TRUE)

##Volcanos-------------------------------------------------------------------------
#quick volcano plot
volcanoplot(fit_bayes, coef=2, highlight = 10, names = rownames(mat),style = "B-statistic")
abline(h = 1.2, v = c(-1.3,1.3), lty = 2)
