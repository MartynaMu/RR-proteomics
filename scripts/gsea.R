# Notes on clusterProfiler------------------------------------------------------
# 1. make sure you have latest version of BiocManager
# 2. update clusterProfiler, enrichplot, org db, dependancies etc
# 3. If still not working, use "remotes::install_github("YuLab-SMU/clusterProfiler"))" 
# to download dev version - last time worked
# 4. Still nothing - download KEGG db and use internal db

# Libs--------------------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)

# geneList prep ----------------------------------------------------------------
# prepare gene List - eg from limma output topTable
# either with adj pvalue cutoff or not

geneList <- deps[,2] #fc
names(geneList) <- as.character(deps[,1]) #gene symbols
geneList <- sort(geneList, decreasing = TRUE) #sort ranks

# In case ud have to switch gene keys
# entrezid <- bitr(as.character(deps[,1]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# 
# geneList <- filter(deps, ID %in% entrezid$SYMBOL)[,2] #fc
# names(geneList) <- as.character(entrezid$ENTREZID) #gene symbols
# geneList <- sort(geneList, decreasing = TRUE) #sort ranks

#write down for external tools, EG. WEBGESTALT
# write.table(geneList, "data/genelist.txt", sep = "\t", quote = FALSE)

# define keytype for genes
keyType = "SYMBOL"

# GO GSEA -------------------------------------------------------------------------
# min/max GSSize means how many genes per term should be used
gsea <- gseGO(geneList = geneList,
              OrgDb = org.Hs.eg.db,
              ont = "BP",
              minGSSize = 3,
              maxGSSize = 100,
              pvalueCutoff = 0.05,
              verbose = TRUE,
              keyType = keyType)

# view results of gsea
# gsea_results <- gsea@result
# upreg <- filter(gsea_results, NES > 0 & p.adjust < 0.05)
# downreg <- filter(gsea_results, NES < 0 & p.adjust < 0.05)

# remove redundant GO terms via GOSemSim methods
simp <- clusterProfiler::simplify(gsea)
simp_results <- simp@result
upreg <- filter(simp_results, NES > 0 & p.adjust < 0.05)
downreg <- filter(simp_results, NES < 0 & p.adjust < 0.05)

# Visualize---------------------------------------------------------------------
# ES - enrichment score - if positive - top of ranked list/upreg, if negative - bottom of ranked list/downreg
# geneSetID - could be a nr or GO code
library(enrichplot)
#gseaplot2(simp, geneSetID = "GO:0022613")

ridgeplot(simp, label_format = 50)+
  geom_vline(xintercept = 0, lty = 2)

## Hierarchical clustering of enriched terms ------------------------------------
# It relies on the pairwise similarities of the enriched terms calculated by the 
# pairwise_termsim() function, which by default using Jaccard’s similarity index (JC). 
# Users can also use semantic similarity values if it is supported (e.g., GO, DO and MeSH).

termism <- pairwise_termsim(gsea)
treeplot(termism)
emapplot(termism)

# 
library(ggplot2)
theme_set(theme_minimal())

upreg %>% 
  arrange(desc(p.adjust)) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  ggplot() +
    geom_col(aes(x = setSize, y = Description, fill = -log10(p.adjust)))+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(axis.text = element_text(size=15))

## Explore expr of terms of interest ------------------------------------------
term = "negative regulation of cell cycle process"
#filtr <- str_split_1(gsea_results$core_enrichment[gsea_results$Description == term], pattern = "/")
filtr <- str_split_1("CDC20/BIRC5/GTSE1/SPC24/SPC25/CCNB1/CDK1/EZH2/CENPF/PPP1R10/NDC80/CDK2/SYF2/BAZ1B/PRPF19/BUB3/NAA10/NPM1/INTS7/H2AX/MAP3K20/RBM14/NOP53/MDC1/TOM1L2/GIGYF2/CDC5L/IK/METTL13/INTS3/WAPL/PSMG2/CHFR/MSH2/BUB1/MAD1L1/CDC73/NUBP1/CDK5RAP2/WDR76/ATRX/TTK/TTI1/ZNF207/ZNF830/NSUN2/TP53/KANK2/TRIAP1
", pattern = "/")

pheatmap::pheatmap(qnorm[qnorm$Gene %in% filtr,1:18],
                   scale = "row",
                   show_rownames = TRUE,
                   clustering_distance_cols = "euclidean",
                   clustering_distance_rows = "correlation",
                   annotation_col = annot_col,
                   annotation_colors = annot_colors,
                   main = term)

# KEGG GSEA-----------------------------------------------------------------------
# works if updated to dev version from github (remotes::install_github("YuLab-SMU/clusterProfiler))
kegg <- clusterProfiler::gseKEGG(geneList = geneList,
                organism = "hsa", 
                keyType = "ncbi-geneid",
                minGSSize = 20,
                maxGSSize = 100,
                pvalueCutoff = 0.05,
                verbose = TRUE)

kegg_result <- kegg@result

ridgeplot(kegg) +
  geom_vline(xintercept = 0, lty = 2)

# Pathways GSEA------------------------------------------------------------------
#doesnt work
pc <- gsePC(geneList = geneList,source = "reactome",keyType = "uniprot")

wp <- gseWP(geneList = geneList,organism = "Homo sapiens")

ridgeplot(wp)
