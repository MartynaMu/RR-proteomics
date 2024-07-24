# Hurdle/MSqrob FC input --------------------------------------------------------
# Assign coefficients
coefs <- seq.int(1,ncol(fc_msqrob),2)
names(coefs) <- fc_msqrob %>% dplyr::select(matches("FC")) %>% colnames() %>% str_replace_all("_FC", "")
coefs

# The input of the following function must be a matrix/dataframe with columns as follows: 
# fc, qval, fc, qval, etc and rownames as gene symbols
##GO GSEA-----------------------------------------------------------------------
compute_gsea <- function(stats, GO = c("BP|CC|MF"), path, coefs, plot = TRUE) {
  ## Case 1. Simplified table results + figures---------
  if (plot == TRUE) {
    # create new paths for results
    path_figs <- paste0(path, "gsego_results/gsego_", GO, "_results/")
    path_data <- paste0(path, "gsego_results/gsego_", GO, "_results/")
    
    if (dir.exists(path_figs) == FALSE) {
      dir.create(path_figs,showWarnings = TRUE,recursive = TRUE)
    }
    
    if (dir.exists(path_data) == FALSE) {
      dir.create(path_data,showWarnings = TRUE,recursive = TRUE)
    }
    # iterate through every comparison [stored as an i, i+1 column in coefs]
    for (i in coefs) {
      # create filenames of results
      filename_data <- paste0("gsego_", GO, "_", names(coefs[coefs==i]), ".tsv")
      filename_fig <- paste0("gsego_", GO, "_top30_", names(coefs[coefs==i]), ".png")
      
      # select the stats for particular comparison, save as a named list
      deps <- stats[c(i,i+1)] %>% filter(get(colnames(stats[i+1])) <= .05)
      geneList <- deps[,1] #fc as a vector
      names(geneList) <- rownames(deps) #gene symbols as names, converted to a list
      geneList <- sort(geneList, decreasing = TRUE) #sort ranks
      
      comp <- str_split_1(names(coefs[coefs==i]), pattern = "_vs_") #prepare exp group names for plot titles
      
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
      if (nrow(gsea_results) > 0) {
        simp <- clusterProfiler::simplify(gsea)
        
        # Visualize and save plots
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
  } else if (plot == FALSE) { #for stat gsea tables only
    ## Case 2. NOT simplified table results ONLY---------
    path_data <- paste0(path, "gsego_results/gsego_", GO, "_results/")
    
    if (dir.exists(path_data) == FALSE) {
      dir.create(path_data,showWarnings = TRUE,recursive = TRUE)
    }
    
    # iterate through every comparison [stored as an i, i+1 column in coefs]
    for (i in coefs) {
      filename_data <- paste0("gsego_", GO, "_", names(coefs[coefs==i]), ".tsv")
      
      # select the stats for particular comparison, save as a named list
      deps <- stats[c(i,i+1)] %>% filter(get(colnames(stats[i+1])) <= .05)
      geneList <- deps[,1] #fc vector
      names(geneList) <- rownames(deps) #gene symbols as names, converted to a list
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
      
      # save gsea stat results
      gsea_results <- gsea@result
      write.table(file = paste0(path_data, filename_data), 
                  x = gsea_results, 
                  quote = FALSE, 
                  sep = "\t",
                  row.names = FALSE)
    }
  }
}

compute_gsea(stats = fc_msqrob, path = "figures/panc1/hurdle_msqrob_fc/", GO = "MF", coefs = coefs, plot = TRUE)
