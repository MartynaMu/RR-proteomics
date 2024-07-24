# x2D_x3Dyoung = x2D - x3D_young, #1
# x2D_x3Dold = x2D - x3D_old, #2
# x2D_PDX2D = x2D - PDX_2D, #3
# x3Dyoung_x3Dold=x3D_young - x3D_old, #4
# x3Dold_PDX3D = x3D_old - PDX_3D, #5
# x3Dyoung_xPDX3D = x3D_young - PDX_3D, #6

# Libs-------------------------------------------------------------------------
library(tidyverse)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

#save coefs for automation
coefs <- seq.int(1,6,1)
names(coefs) <- c("2dv3dy", "2dv3do", "2dvPDX2d", "3dyv3do", "3dovPDX3d", "3dyvPDX3d")

#Automatized--------------------------------------------------------------------
##GO GSEA-----------------------------------------------------------------------
compute_gsea <- function(fit, GO = c("BP|CC|MF"), coefs, plot = TRUE) {
  ## Case 1. Simplified table results + figures---------
  if (plot == TRUE) {
    path_figs <- paste0("figures/gsego_results/gsego_", GO, "_results/")
    path_data <- paste0("data/gsego_results/gsego_", GO, "_results/")
      
    if (dir.exists(path_figs) == FALSE) {
      dir.create(path_figs,showWarnings = TRUE,recursive = TRUE)
    }
    
    if (dir.exists(path_data) == FALSE) {
      dir.create(path_data,showWarnings = TRUE,recursive = TRUE)
    }
    
    for (i in coefs) {
      filename_data <- paste0("gsego_", GO, "_", names(coefs[i]), ".tsv")
      filename_fig <- paste0("gsego_", GO, "_top30_", names(coefs[i]), ".png")
      
      deps <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(qnorm))
  
      geneList <- deps[,2] #fc
      names(geneList) <- as.character(deps[,1]) #gene symbols
      geneList <- sort(geneList, decreasing = TRUE) #sort ranks
      
      comp <- str_split_1(names(coefs[i]), pattern = "v") #prepare exp group names for plot titles
      
      # define keytype for genes
      keyType = "SYMBOL"
      
      # GO GSEA 
      # min/max GSSize means how many genes per term should be used
      gsea <- gseGO(geneList = geneList,
                    OrgDb = org.Hs.eg.db,
                    ont = GO,
                    minGSSize = 15,
                    maxGSSize = 100,
                    pvalueCutoff = 0.05,
                    verbose = TRUE,
                    keyType = keyType)
    
      
      # write down unsimplified results
      gsea_results <- gsea@result
      write.table(file = paste0(path_data, filename_data), 
                  x = gsea_results, 
                  quote = FALSE, 
                  sep = "\t",
                  row.names = FALSE)
      
      # remove redundant GO terms via GOSemSim methods for plotting
      simp <- clusterProfiler::simplify(gsea)
      
      # Visualize
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
  } else if (plot == FALSE) { #for tables only
    ## Case 2. NOT simplified table results ONLY---------
      path_data <- paste0("data/gsego_results/gsego_", GO, "_results/")
      
      if (dir.exists(path_data) == FALSE) {
        dir.create(path_data,showWarnings = TRUE,recursive = TRUE)
      }
      
      for (i in coefs) {
        filename_data <- paste0("gsego_", GO, "_", names(coefs[i]), ".tsv")
        
        deps <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(qnorm))
        
        geneList <- deps[,2] #fc
        names(geneList) <- as.character(deps[,1]) #gene symbols
        geneList <- sort(geneList, decreasing = TRUE) #sort ranks
        
        # define keytype for genes
        keyType = "SYMBOL"
        
        # GO GSEA 
        # min/max GSSize means how many genes per term should be used
        gsea <- gseGO(geneList = geneList,
                      OrgDb = org.Hs.eg.db,
                      ont = GO,
                      minGSSize = 15,
                      maxGSSize = 100,
                      pvalueCutoff = 0.05,
                      verbose = TRUE,
                      keyType = keyType)
        
        
        gsea_results <- gsea@result
        write.table(file = paste0(path_data, filename_data), 
                    x = gsea_results, 
                    quote = FALSE, 
                    sep = "\t",
                    row.names = FALSE)
    }
  }
}

compute_gsea(fit = fit_bayes, GO = "MF", coefs = coefs, plot = TRUE)


## KEGG GSEA-------------------------------------------------------------------
compute_keggsea <- function(fit, coefs) {
  path_figs <- paste0("figures/gsekegg_results/")
  path_data <- paste0("data/gsekegg_results/")
  if (dir.exists(path_figs) == FALSE) {
    dir.create(path_figs,showWarnings = TRUE,recursive = TRUE)
  }
  
  if (dir.exists(path_data) == FALSE) {
    dir.create(path_data,showWarnings = TRUE,recursive = TRUE)
  }
  
  for (i in coefs) {
    filename_data <- paste0("gsekegg_", names(coefs[i]), ".tsv")
    filename_fig <- paste0("gsekegg_", "top30_", names(coefs[i]), ".png")
    
    deps <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(qnorm))
    
    entrezid <- bitr(as.character(deps[,1]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    geneList <- filter(deps, ID %in% entrezid$SYMBOL)[,2] #fc
    names(geneList) <- as.character(entrezid$ENTREZID) #gene symbols
    geneList <- sort(geneList, decreasing = TRUE) #sort ranks
    
    comp <- str_split_1(names(coefs[i]), pattern = "v") #prepare exp group names for plot titles
    
    kegg <- gseKEGG(geneList = geneList,                
                    organism = "hsa", 
                    keyType = "ncbi-geneid",
                    minGSSize = 20,
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

compute_keggsea(fit=fit_bayes, coefs = coefs)

# Manual --------------------------------------------------------------------------
deps <- topTable(fit_bayes, coef = 2, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(qnorm))

# entrezid <- bitr(as.character(deps[,1]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# geneList <- filter(deps, ID %in% entrezid$SYMBOL)[,2] #fc
# names(geneList) <- as.character(entrezid$ENTREZID) #gene symbols
# geneList <- sort(geneList, decreasing = TRUE) #sort ranks

geneList <- deps[,2] #fc
names(geneList) <- as.character(deps[,1]) #gene symbols
geneList <- sort(geneList, decreasing = TRUE) #sort ranks

# define keytype for genes
keyType = "ENTREZID"
keyType = "SYMBOL"

# GO GSEA 
# min/max GSSize means how many genes per term should be used
gsea <- gseGO(geneList = geneList,
              OrgDb = org.Hs.eg.db,
              ont = "BP",
              minGSSize = 15,
              maxGSSize = 100,
              pvalueCutoff = 0.05,
              verbose = TRUE,
              keyType = keyType)

kegg <- gseKEGG(geneList = geneList,                
                organism = "hsa", 
                keyType = "ncbi-geneid",
                minGSSize = 20,
                maxGSSize = 100,
                pvalueCutoff = 0.05,
                verbose = TRUE)
kegg_result <- kegg@result

# remove redundant GO terms via GOSemSim methods
simp <- clusterProfiler::simplify(gsea)
simp_results <- simp@result
write.table(file = paste0("data/gsego_results/gsego_BP_results/unfiltered/gsego_bp_",
                          names(coefs[1]),".tsv"), x = simp_results, quote = FALSE, sep = "\t")

# Visualize
ridgeplot(kegg, label_format = 60)+
  geom_vline(xintercept = 0, lty = 2)+
  labs(x = "log2FC")

ggsave(filename = paste0("gsego_bp_top30_", names(coefs[1]), ".png"),
       path = "figures/gsego_results/gsego_BP_results/unfiltered",
       device = "png",
       width = 800,
       height = 780,
       units = "px",
       dpi = 100)

