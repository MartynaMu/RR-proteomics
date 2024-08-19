# Libs-------------------------------------------------------------------------
library(tidyverse)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Intensity mat FC input ------------------------------------------------------------
coefs <- seq.int(1,18,1)
names(coefs) <- colnames(fit_bayes$contrasts)

## GO GSEA-----------------------------------------------------------------------
compute_gsea <- function(fit, path, GO = c("BP|CC|MF"), coefs) {
  path_figs <- paste0(path, "gsego_results/gsego_", GO, "_results/")
  path_data <- paste0(path, "gsego_results/gsego_", GO, "_results/")
  
  if (dir.exists(path_figs) == FALSE) {
    dir.create(path_figs,showWarnings = TRUE,recursive = TRUE)
  }
  
  if (dir.exists(path_data) == FALSE) {
    dir.create(path_data,showWarnings = TRUE,recursive = TRUE)
  }
  
  for (i in coefs) {
    filename_data <- paste0("gsego_", GO, "_", names(coefs[i]), ".tsv")
    filename_fig <- paste0("gsego_", GO, "_top30_", names(coefs[i]), ".png")
    
    deps <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(mat))
    
    geneList <- deps[,2] #fc
    names(geneList) <- as.character(deps[,1]) #gene symbols
    geneList <- sort(geneList, decreasing = TRUE) #sort ranks
    
    comp <- str_split_1(names(coefs[i]), pattern = "_") #prepare exp group names for plot titles
    
    # define keytype for genes
    keyType = "SYMBOL"
    
    # GO GSEA 
    # min/max GSSize means how many genes per term should be used
    gsea <- gseGO(geneList = geneList,
                  OrgDb = org.Hs.eg.db,
                  ont = GO,
                  minGSSize = 5,
                  maxGSSize = 100,
                  pvalueCutoff = 0.05,
                  verbose = TRUE,
                  keyType = keyType)
    
    
    # write down unsimplified results
    gsea_results <- gsea@result
    write.table(file = paste0(path_data, "unsimplified.", filename_data), 
                x = gsea_results, 
                quote = FALSE, 
                sep = "\t",
                row.names = FALSE)
    
    # remove redundant GO terms via GOSemSim methods for plotting
    simp <- clusterProfiler::simplify(gsea)
    
    # write down simplified results
    write.table(file = paste0(path_data, "simplified.", filename_data), 
                x = simp@result, 
                quote = FALSE, 
                sep = "\t",
                row.names = FALSE)
    
    # Visualize
    if (nrow(gsea_results) > 0) {
      ridgeplot(simp, label_format = 60)+
        geom_vline(xintercept = 0, lty = 2)+
        labs(x = sprintf("log2 FC (%s/%s)", comp[2], comp[1]))
      
      ggsave(filename = filename_fig,
             device = "png",
             path = path_figs,
             width = 800,
             height = 780,
             units = "px",
             dpi = 100)
    }
  }
}

compute_gsea(fit = fit_bayes, path = "figures/allruns/final_quant/", GO = "MF", coefs = coefs)

## KEGG GSEA-------------------------------------------------------------------
compute_keggsea <- function(fit, path, coefs) {
  path_figs <- paste0(path, "gsekegg_results")
  path_data <- paste0(path, "gsekegg_results")
  if (dir.exists(path_figs) == FALSE) {
    dir.create(path_figs,showWarnings = TRUE,recursive = TRUE)
  }
  
  if (dir.exists(path_data) == FALSE) {
    dir.create(path_data,showWarnings = TRUE,recursive = TRUE)
  }
  
  for (i in coefs) {
    filename_data <- paste0("gsekegg_", names(coefs[i]), ".tsv")
    filename_fig <- paste0("gsekegg_", "top30_", names(coefs[i]), ".png")
    
    deps <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(mat))
    
    entrezid <- bitr(as.character(deps[,1]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    geneList <- filter(deps, ID %in% entrezid$SYMBOL)[,2] #fc
    names(geneList) <- as.character(entrezid$ENTREZID) #gene symbols
    geneList <- sort(geneList, decreasing = TRUE) #sort ranks
    
    comp <- str_split_1(names(coefs[i]), pattern = "_") #prepare exp group names for plot titles
    
    kegg <- gseKEGG(geneList = geneList,                
                    organism = "hsa", 
                    keyType = "ncbi-geneid",
                    minGSSize = 5,
                    maxGSSize = 100,
                    pvalueCutoff = 0.05,
                    verbose = TRUE)
    
    kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    
    kegg_results <- kegg@result
    write.table(file = paste0(path_data, filename_data), 
                x = kegg_results, 
                quote = FALSE, 
                sep = "\t",
                row.names = FALSE)
    
    # Visualize
    if (nrow(kegg_results) > 0) {
      ridgeplot(kegg, label_format = 60)+
        geom_vline(xintercept = 0, lty = 2)+
        labs(x = sprintf("log2 FC (%s/%s)", comp[2], comp[1]))
      
      ggsave(filename = filename_fig,
             device = "png",
             path = path_figs,
             width = 800,
             height = 780,
             units = "px",
             dpi = 100)
    }
  }
}

compute_keggsea(fit = fit_bayes, path = "figures/allruns/final_quant/", coefs = coefs)

## WikiPathways ------------------------------------------------
compute_WPgsea <- function(fit, path, coefs) {
  path_figs <- paste0(path, "gseWP_results")
  path_data <- paste0(path, "gseWP_results")
  if (dir.exists(path_figs) == FALSE) {
    dir.create(path_figs,showWarnings = TRUE,recursive = TRUE)
  }
  
  if (dir.exists(path_data) == FALSE) {
    dir.create(path_data,showWarnings = TRUE,recursive = TRUE)
  }
  
  for (i in coefs) {
    filename_data <- paste0("gseWP_", names(coefs[i]), ".tsv")
    filename_fig <- paste0("gseWP_", "top30_", names(coefs[i]), ".png")
    
    deps <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(mat))
    
    entrezid <- bitr(as.character(deps[,1]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    geneList <- filter(deps, ID %in% entrezid$SYMBOL)[,2] #fc
    names(geneList) <- as.character(entrezid$ENTREZID) #gene symbols
    geneList <- sort(geneList, decreasing = TRUE) #sort ranks
    
    comp <- str_split_1(names(coefs[i]), pattern = "_") #prepare exp group names for plot titles
    
    wp <- gseWP(geneList = geneList,                
                 organism = "Homo sapiens")
    
    wp <- setReadable(wp, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    
    wp_results <- wp@result
    write.table(file = paste0(path_data, filename_data), 
                x = wp_results, 
                quote = FALSE, 
                sep = "\t",
                row.names = FALSE)
    
    # Visualize
    if (nrow(wp_results) > 0) {
      ridgeplot(wp, label_format = 60)+
        geom_vline(xintercept = 0, lty = 2)+
        labs(x = sprintf("log2 FC (%s/%s)", comp[2], comp[1]))
      
      ggsave(filename = filename_fig,
             device = "png",
             path = path_figs,
             width = 800,
             height = 780,
             units = "px",
             dpi = 100)
    }
  }
}

compute_WPgsea(fit = fit_bayes, path = "figures/allruns/final_quant/", coefs = coefs)
