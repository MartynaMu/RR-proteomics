cell.line <- c("miapaca", "panc1", "cfpac")

gene_occ <- data.frame(matrix(ncol = 6))
for (l in cell.line) {
  library(tidyverse)
  fc_msqrob <- read.delim(sprintf("data/%s/a_Hurdle_Msqrob.tsv",cell.line))
  # fc_msqrob <- read.delim("data/panc1/a_msqrob_result_Ridge_Aggregate.tsv")
  
  # Preprocessing-------------------------------------------------------------------
  #tidy Gene symbols
  fc_msqrob$Genes <- str_replace_all(fc_msqrob$Genes, "_HUMAN", "")
  colnames(fc_msqrob) <- str_replace_all(colnames(fc_msqrob), "norm_", "")
  
  prot_info <- fc_msqrob %>% dplyr::select(matches("Gene|Protein|Peptide")) # protein info
  int_mat <- fc_msqrob %>% dplyr::select(!matches("FC|qval|_OR|_fisher|_pval")) # intensity matrix
  colnames(int_mat)
  
  #first check filtered on normal qval, then fisher_qval
  fc_msqrob <- fc_msqrob %>% dplyr::select(matches("FC|(?<!fisher)_qval|Genes", perl = TRUE)) 
  #fc_msqrob <- fc_msqrob %>% dplyr::select(matches("FC|(?<=fisher)_qval|Genes", perl = TRUE)) 
  fc_msqrob <- column_to_rownames(fc_msqrob, "Genes")
  colnames(fc_msqrob)
  
  fc_msqrob <- fc_msqrob[,order(colnames(fc_msqrob))]
  colnames(fc_msqrob)
  
  # get rid of 0
  int_mat[int_mat == 0] <- NA
  
  # normalize quantiles
  int_mat <- column_to_rownames(int_mat, "Genes")
  int_mat <- int_mat %>% dplyr::select(!matches("Protein|Peptide"))
  colnames(int_mat)
  qnorm <- limma::normalizeQuantiles(drop_na(int_mat))
  
  # DE limma --------------------------------------------------------------------
  # limma instead of hurdle/msqrob ttest
  
  library(tidyverse)
  library(limma)
  
  # In-depth comparisons-----------------------------------------------------------
  ##Design table creation-------------------------------------------------------------
  #where first digit in rep represents the column it will be filled with, second is nr of replicates
  
  mat <- as.matrix(qnorm[1:ncol(qnorm)])
  V1 <- str_count(colnames(mat), pattern = "2D(?!_XENO)") %>% sum()
  V2 <- str_count(colnames(mat), pattern = "3D_YOUNG") %>% sum()
  V3 <- str_count(colnames(mat), pattern = "3D_OLD") %>% sum()
  V4 <- str_count(colnames(mat), pattern = "2D_XENO") %>% sum()
  V5 <- str_count(colnames(mat), pattern = "3D_XENO") %>% sum()
  
  design <- model.matrix(~ 0+factor(c(rep(1,V1),
                                      rep(2,V2),
                                      rep(3,V3),
                                      rep(4,V4),
                                      rep(5,V5))))
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
  
  
  # GSEA ------------------------------------------------------------
  ## Libs
  library(tidyverse)
  library(limma)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  
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
  
  test <- pairwise_termsim(gsea)
  test <- simplify(test, cutoff=0.7,by="p.adjust", select_fun = min)
  
  
  # count terms that every gene was taken account into
  gene_occ <- data.frame(matrix(ncol=2))
  colnames(gene_occ) <- c("Gene", "TermCount")
  for (i in rownames(qnorm)) {
    temp <- test@result$core_enrichment |> str_count(pattern = i) %>% sum()
    gene_occ <- add_row(gene_occ, Gene = i, TermCount = temp)
  }
  gene_occ <- gene_occ %>% drop_na() %>% arrange(desc(TermCount))
  head(gene_occ)
  assign(paste("gene_occ",cell.line,names(coefs[1]),sep="."), gene_occ)
}

test <- full_join(gene_occ.cfpac.2dv3dy, gene_occ.miapaca.2dv3dy, by="Gene") |> full_join(gene_occ.panc1.2dv3dy, by="Gene")
  
test |> pivot_longer(2:4, names_to = "CellLine", values_to = "TermCount") |> 
  filter(TermCount > 7) %>%
  ggplot(aes(x=Gene,y=TermCount))+
  geom_col()

## Reactome ------------------------------------------------------------------

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




# Visualize -------------------------------------------------------------------

# Visualize
ridgeplot(simp, label_format = 60)+
  geom_vline(xintercept = 0, lty = 2)+
  labs(x = sprintf("log2 FC (%s/%s)", comp[2], comp[1]))

cnetplot(filter(test, NES > 0), 
         foldChange = geneList, 
         categorySize = "geneNum",
         node_label = "gene",
         layout = "graphopt",
         cex_label_gene = .5, 
         showCategory = 108)

ggsave("plot.pdf", width =1000, height =1000, units = "mm")

emapplot(test, cex_label_category = .7, repel= TRUE, showCategory = 138)

treeplot(test, showCategory = nrow(simp_results), nCluster = 15)

heatplot(test, foldChange = geneList, pvalue = "pvalue", showCategory = 138)