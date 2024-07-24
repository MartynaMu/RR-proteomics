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

# Automatically--------------------------------------------------------------------
enrich_go <- function(fit, coefs, GO="BP|CC|MF", simplify = TRUE) {
  # simplify=TRUE -------------------------
  if (simplify == TRUE) {
    path_figs <- paste0("figures/enrichgo_results/enrich_", GO, "_results/unfiltered/")
    path_data <- paste0("data/enrichgo_results/enrich_", GO, "_results/unfiltered/simplified/")
    
    if (dir.exists(path_figs) == FALSE) {
      dir.create(path_figs,showWarnings = TRUE,recursive = TRUE)
    }
    
    if (dir.exists(path_data) == FALSE) {
      dir.create(path_data,showWarnings = TRUE,recursive = TRUE)
    }
    
    for (i in coefs) {
      filename_data <- paste0("enrich_", GO, "_sim_", names(coefs[i]))
      filename_fig <- paste0("enrich_", GO, "_top30_sim_", names(coefs[i]))
      
      deps <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(qnorm))
      
      geneList <- deps[,2] #fc
      names(geneList) <- as.character(deps[,1]) #gene symbols
      geneList <- sort(geneList, decreasing = TRUE) #sort ranks
      
      comp <- str_split_1(names(coefs[i]), pattern = "v") #prepare exp group names for plot titles
      
      # define keytype for genes
      keyType = "SYMBOL"
      
      #extract ids
      upreg <- deps$ID[deps$logFC >= 1.3 & deps$adj.P.Val <= 0.05]
      downreg <- deps$ID[deps$logFC <= -1.3 & deps$adj.P.Val <= 0.05]
      message(paste0(names(coefs[i]),": ", "Nr of upregulated genes:", length(upreg)))
      message(paste0(names(coefs[i]),": ", "Nr of downregulated genes:", length(downreg)))
      
      #enrich upreg
      go_enr <- enrichGO(gene = upreg,
                         universe = names(geneList),
                         OrgDb = org.Hs.eg.db,
                         keyType = "SYMBOL",
                         ont = GO,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH")
      go_enr_results <- go_enr@result
      write.table(file = paste0(path_data, filename_data, "_upreg.tsv"), 
                  x = go_enr_results, 
                  quote = FALSE, 
                  sep = "\t",
                  dec = ",",
                  row.names = FALSE)
      #Visualize upreg
      if(min(go_enr_results$p.adjust) <= 0.05) {
        go_enr_sim <- simplify(go_enr)
        barplot(go_enr_sim, showCategory = 30, label_format = 60)+
          labs(title = sprintf("Upregulated genes in %s",comp[2]))
        ggsave(filename = paste0(filename_fig,"_upreg.png"),
               device = "png",
               path = path_figs,
               width = 800,
               height = 780,
               units = "px",
               dpi = 100)
      } else {
        message(paste0("No enrichment for upregulated genes in ", names(coefs[i])))
      }
      #enrich downreg
      go_enr <- enrichGO(gene = downreg,
                         universe = names(geneList),
                         OrgDb = org.Hs.eg.db,
                         keyType = "SYMBOL",
                         ont = GO,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH")
      go_enr_results <- go_enr@result
      write.table(file = paste0(path_data, filename_data, "_downreg.tsv"), 
                  x = go_enr_results, 
                  quote = FALSE, 
                  sep = "\t",
                  dec = ",",
                  row.names = FALSE)
      #Visualize downreg
      if(min(go_enr_results$p.adjust) <= 0.05) {
        go_enr_sim <- simplify(go_enr)
        barplot(go_enr_sim, showCategory = 30, label_format = 60)+
          labs(title = sprintf("Downregulated genes in %s",comp[2]))
        ggsave(filename = paste0(filename_fig,"_downreg.png"),
               device = "png",
               path = path_figs,
               width = 800,
               height = 780,
               units = "px",
               dpi = 100)
      } else {
        message(paste0("No enrichment for downregulated genes in ", names(coefs[i])))
      }
      # simplify=FALSE ----------------------------------------------
    } 
  } else if (simplify == FALSE) {
    path_data <- paste0("data/enrichgo_results/enrich_", GO, "_results/unfiltered/not_simplified/")
    
    if (dir.exists(path_figs) == FALSE) {
      dir.create(path_figs,showWarnings = TRUE,recursive = TRUE)
    }
    
    if (dir.exists(path_data) == FALSE) {
      dir.create(path_data,showWarnings = TRUE,recursive = TRUE)
    }
    
    for (i in coefs) {
      filename_data <- paste0("enrich_", GO, "_notsim_", names(coefs[i]))
      
      deps <- topTable(fit_bayes, coef = i, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(qnorm))
      
      geneList <- deps[,2] #fc
      names(geneList) <- as.character(deps[,1]) #gene symbols
      geneList <- sort(geneList, decreasing = TRUE) #sort ranks
      
      # define keytype for genes
      keyType = "SYMBOL"
      
      #extract ids
      upreg <- deps$ID[deps$logFC >= 1.3 & deps$adj.P.Val <= 0.05]
      downreg <- deps$ID[deps$logFC <= -1.3 & deps$adj.P.Val <= 0.05]
      message(paste0(names(coefs[i]),": ", "Nr of upregulated genes:", length(upreg)))
      message(paste0(names(coefs[i]),": ", "Nr of downregulated genes:", length(downreg)))
      
      #enrich upreg
      go_enr <- enrichGO(gene = upreg,
                         universe = names(geneList),
                         OrgDb = org.Hs.eg.db,
                         keyType = "SYMBOL",
                         ont = GO,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH")
      go_enr_results <- go_enr@result
      write.table(file = paste0(path_data, filename_data, "_notsim_upreg.tsv"), 
                  x = go_enr_results, 
                  quote = FALSE, 
                  sep = "\t",
                  dec = ",",
                  row.names = FALSE)

      if(min(go_enr_results$p.adjust) >= 0.05) {
        message(paste0("No enrichment for upregulated genes in ", names(coefs[i])))
      } 
      #enrich downreg
      go_enr <- enrichGO(gene = downreg,
                         universe = names(geneList),
                         OrgDb = org.Hs.eg.db,
                         keyType = "SYMBOL",
                         ont = GO,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH")
      go_enr_results <- go_enr@result
      write.table(file = paste0(path_data, filename_data, "_notsim_downreg.tsv"), 
                  x = go_enr_results, 
                  quote = FALSE, 
                  sep = "\t",
                  dec = ",",
                  row.names = FALSE)
      if(min(go_enr_results$p.adjust) >= 0.05) {
        message(paste0("No enrichment for downregulated genes in ", names(coefs[i])))
      }
    }
  }  
}
  
enrich_go(fit_bayes, coefs=coefs, GO="MF", simplify = TRUE)

# Manually---------------------------------------------------------------------
deps <- topTable(fit_bayes, coef = 2, adjust="BH", genelist = rownames(mat), resort.by = "logFC", number = nrow(qnorm))

geneList <- deps[,2] #fc
names(geneList) <- as.character(deps[,1]) #gene symbols
geneList <- sort(geneList, decreasing = TRUE) #sort ranks

# define keytype for genes
keyType = "SYMBOL"

#extract ids
upreg <- deps$ID[deps$logFC >= 1.3 & deps$adj.P.Val <= 0.05]
downreg <- deps$ID[deps$logFC <= -1.3 & deps$adj.P.Val <= 0.05]

#enrich upreg
go_enr <- enrichGO(gene = downreg,
                   universe = names(geneList),
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",
                   pvalueCutoff = 0.05)
go_enr_results <- go_enr@result
write.table(file = paste0("data/enrichgo_results/enrich_BP_results/unfiltered/enrich_BP_", names(coefs[2]), "_upreg.tsv"), 
            x = go_enr_results, 
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE)
#Visualize upreg
barplot(go_enr, showCategory = 30, label_format = 60)
ggsave(filename = paste0("enrich_BP_top30_", names(coefs[2]),"_upreg.png"),
       device = "png",
       path = path_figs,
       width = 800,
       height = 780,
       units = "px",
       dpi = 100)

#enrich downreg
go_enr <- enrichGO(gene = downreg,
                   universe = names(geneList),
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",
                   qvalueCutoff = 0.05,
                   pAdjustMethod = "BH")
go_enr_results <- go_enr@result
write.table(file = paste0(path_data, filename_data, "_downreg.tsv"), 
            x = go_enr_results, 
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE)
#Visualize downreg
barplot(go_enr, showCategory = 30, label_format = 60)
ggsave(filename = paste0(filename_fig,"_downreg.png"),
       device = "png",
       path = path_figs,
       width = 800,
       height = 780,
       units = "px",
       dpi = 100)