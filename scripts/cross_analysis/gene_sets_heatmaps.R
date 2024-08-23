library(tidyverse)
library(clusterProfiler)

# Retrieve gene symbols from overlapping terms -------------------------------
# 1. prepare a list of pathways enriched in at least 2 cell lines
overlap.wp <- unlist(full.overlap.terms) %>% unique()

# 2. Find each pathway in gsea objects and retrieve info from the pathway
term.desc <- data.frame(matrix(ncol=3,nrow=1))
colnames(term.desc) <- c("ID", "Description", "core_enrichment")
for (i in wp_list) {
  term.desc <- bind_rows(term.desc, i@result %>% filter(ID %in% overlap.wp) %>% select(ID, Description, core_enrichment))
}
term.desc <- term.desc %>% group_by(ID, Description) %>% summarise(Genes = str_flatten(pull(across(core_enrichment))))
term.desc <- term.desc %>% rowwise() %>% mutate(Genes = str_split(Genes, pattern = "/"))
term.desc <- drop_na(term.desc)

# Prepare a matrix of FC and p-values from each comparison in each cell line ------------
args <- c(number = nrow(mat), resort.by = "logFC")

temp <- data.frame("Genes" = rownames(topTable(fit_bayes, number=nrow(mat))))
for (i in seq.int(1,6,1)) {
  temp <- left_join(temp,(topTable(fit_bayes,coef=i,substitute(args))) |> select(logFC) %>% rownames_to_column(var = "Genes"), by="Genes")
  temp <- left_join(temp,(topTable(fit_bayes,coef=i+6,substitute(args))) |> select(logFC) %>% rownames_to_column(var = "Genes"), by="Genes")
  temp <- left_join(temp,(topTable(fit_bayes,coef=i+12,substitute(args))) |> select(logFC) %>% rownames_to_column(var = "Genes"), by="Genes")
  colnames(temp)[(ncol(temp)-2):ncol(temp)] <- str_c(cell.line,curr.comp[i],"logFC", sep = ".")
}
assign("FC", temp)
FC <- column_to_rownames(FC, var="Genes")

# Prepare a matrix of means in each comparison in each cell line -------------------------------
means.col <- colnames(test) %>% str_sub(end=-3L) %>% unique()
gr.means <- data.frame(matrix(ncol=1,nrow=nrow(test)))
for (i in means.col) {
 gr.means <- bind_cols(gr.means,test %>% rowwise() %>% 
          transmute(mean(c_across(matches(paste0(i,"_[1234]"))), na.rm=TRUE)))
          
}
gr.means <- gr.means[-1]
colnames(gr.means) <- means.col
rownames(gr.means) <- rownames(test)
# order colnames group wise
order <- colnames(gr.means) %>% str_replace_all(pattern = "MIAPACA|PANC|CFPAC", replacement = "") %>% unique()
order[1] <- str_c(order[1],"$")

temp <- c()
for (i in order) {
  temp <- c(temp, str_which(colnames(gr.means), pattern = i))
}

gr.means <- relocate(gr.means, all_of(temp))

# 5. Retrieve FC of the genes from each comparison in each cell line
FC_wp <- FC %>% filter(rownames(FC) %in% term.genes.vec)

# 6. Look for a specific term genes
term <- term.genes[["WP2363"]]


# 7. Draw heatmap
# HM loop for all terms --------------------------------------------------------------

terms.hm <- function(FC.matrix = NULL, mean.matrix = NULL, term.descr, save = TRUE, path = NULL) {
  
  if (!is.null(FC.matrix)) {
    matrix <- FC.matrix
    plot.desc <- "Log2FC in each comparison"
    file.desc <- ".hm.fc.png"
    path2 <- "fc/"
  } else if (!is.null(mean.matrix)) {
    matrix <- mean.matrix
    plot.desc <- "Mean group intensity"
    file.desc <- ".hm.mean.png"
    path2 <- "means/"
  }
  terms <- term.descr$ID
  for (i in terms) {
    term.desc <- term.descr %>% filter(ID == i) %>% pull(Description)
    term.genes <- term.descr %>% filter(ID == i) %>% pull(Genes)
    names(term.genes) <- i
    matrix.to.use <- matrix %>% filter(rownames(matrix) %in% term.genes[[i]])
    if (length(term.genes[[i]]) > 20) {
      p <- pheatmap::pheatmap(matrix.to.use, 
                              cluster_cols = FALSE,
                              gaps_col = c(3,6,9,12,15),
                              cluster_rows = TRUE,
                              angle_col = 45,
                              border_color = "black",
                              main = sprintf("%s - %s (%s)", plot.desc, term.desc, i),
                              col=colorRampPalette(c("cornflowerblue","white","red"))(100))
      if (save==TRUE) {
        ggsave(filename = paste0(i, file.desc),
               plot = p,
               path = paste0(path ,path2),
               device = "png",
               height = 500,
               width = 1100,
               units = "px",
               dpi = 100)
      }
    } else if (length(term.genes[[i]]) <= 20) {
      p <- pheatmap::pheatmap(matrix.to.use, 
                              cluster_cols = FALSE,
                              gaps_col = c(3,6,9,12,15),
                              cluster_rows = TRUE,
                              angle_col = 45,
                              display_numbers = TRUE,
                              fontsize_number = 10,
                              number_color = "black",
                              border_color = "black",
                              main = sprintf("%s - %s (%s)", plot.desc, term.desc, i),
                              col=colorRampPalette(c("cornflowerblue","white","red"))(100))
      if (save==TRUE) {
        ggsave(filename = paste0(i, file.desc),
               plot = p,
               path = paste0(path, path2),
               device = "png",
               height = 500,
               width = 1100,
               units = "px",
               dpi = 100)
      }
    }
  }
}

terms.hm(FC.matrix = FC_wp, term.descr = term.desc, save = TRUE, path = "figures/allruns/final_quant/term_overlap/")

# KEGG panc cancer gene set tested ---------------------------------------------
# Parse GMT file for KEGG pancreatic cancer hallmark genes

#temp <- read.gmt(gmtfile = "data/KEGG_PANCREATIC_CANCER.v2024.1.Hs.gmt") # this one doesn't work for now
temp <- read.gmt.wp(gmtfile = "data/KEGG_PANCREATIC_CANCER.v2024.1.Hs.gmt")

