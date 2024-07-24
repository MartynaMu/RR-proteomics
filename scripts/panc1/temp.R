# Libs-------------------------------------------------------------------------
library(tidyverse)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Intensity mat FC input ------------------------------------------------------------
coefs <- seq.int(1,6,1)
names(coefs) <- c("2dv3dy", "2dv3do", "2dvPDX2d", "3dyv3do", "3dovPDX3d", "3dyvPDX3d")

      
deps <- topTable(fit_bayes, coef = 1, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(int_mat))

geneList <- deps[,2] #fc
names(geneList) <- as.character(deps[,1]) #gene symbols
geneList <- sort(geneList, decreasing = TRUE) #sort ranks

comp <- str_split_1(names(coefs[1]), pattern = "v") #prepare exp group names for plot titles

# define keytype for genes
keyType = "SYMBOL"

# GO GSEA 
# min/max GSSize means how many genes per term should be used
gsea <- gseGO(geneList = geneList,
              OrgDb = org.Hs.eg.db,
              ont = "BP",
              minGSSize = 1,
              maxGSSize = 100,
              pvalueCutoff = 0.05,
              verbose = TRUE,
              keyType = keyType)


# write down unsimplified results
gsea_results <- gsea@result

# remove redundant GO terms via GOSemSim methods for plotting
simp <- clusterProfiler::simplify(gsea)
simp <- clusterProfiler::simplify(simp)
simp_results <- simp@result

# Visualize
ridgeplot(simp, label_format = 60)+
    geom_vline(xintercept = 0, lty = 2)+
    labs(x = sprintf("log2 FC (%s/%s)", comp[2], comp[1]))
test <- pairwise_termsim(gsea)
test <- simplify(test, cutoff=0.7,by="p.adjust", select_fun = min)

cnetplot(test, foldChange = geneList, categorySize = "pvalue", cex_label_gene = .5, showCategory = 10)

emapplot(test, cex_label_category = .7, repel= TRUE, showCategory = 138)

treeplot(test, showCategory = nrow(simp_results), nCluster = 15)

heatplot(test, foldChange = geneList, pvalue = "pvalue", showCategory = 138)

library(ReactomePA)


## KEGG GSEA-------------------------------------------------------------------
deps <- topTable(fit_bayes, coef = 1, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(int_mat))

entrezid <- bitr(as.character(deps[,1]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
geneList <- filter(deps, ID %in% entrezid$SYMBOL)[,2] #fc
names(geneList) <- as.character(entrezid$ENTREZID) #gene symbols
geneList <- sort(geneList, decreasing = TRUE) #sort ranks

comp <- str_split_1(names(coefs[1]), pattern = "v") #prepare exp group names for plot titles

kegg <- gseKEGG(geneList = geneList,                
                organism = "hsa", 
                keyType = "ncbi-geneid",
                minGSSize = 20,
                maxGSSize = 100,
                pvalueCutoff = 0.05,
                verbose = TRUE)

kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

kegg_results <- kegg@result
    
# Visualize
ridgeplot(kegg, label_format = 60)+
    geom_vline(xintercept = 0, lty = 2)+
    labs(x = sprintf("log2 FC (%s/%s)", comp[2], comp[1]))
  
