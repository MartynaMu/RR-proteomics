library(tidyverse)
library(clusterProfiler)

# Retrieve gene symbols from overlapping terms -------------------------------
# 1. prepare a list of pathways enriched in at least 2 cell lines
overlap.wp <- unlist(full.overlap.terms) %>% unique()
#overlap.kegg <- unlist(full.overlap.terms) %>% unique()
#overlap.bp <- unlist(full.overlap.terms) %>% unique()
#overlap.cc <- unlist(full.overlap.terms) %>% unique()
#overlap.mf <- unlist(full.overlap.terms) %>% unique()

# 2. Find each pathway in gsea objects and retrieve info from the pathway
term.desc <- data.frame(matrix(ncol=3,nrow=1))
colnames(term.desc) <- c("ID", "Description", "core_enrichment")
for (i in wp_list) {
  term.desc <- bind_rows(term.desc, i@result %>% filter(ID %in% overlap.wp) %>% dplyr::select(ID, Description, core_enrichment))
}
term.desc <- term.desc %>% group_by(ID, Description) %>% summarise(Genes = str_flatten(pull(across(core_enrichment)),collapse = "/"))
term.desc <- term.desc %>% rowwise() %>% mutate(Genes = str_split(Genes, pattern = "/") |> map(unique))
term.desc <- drop_na(term.desc)

#run in case of GO gsea
term.desc <- term.desc %>% mutate(ID = str_sub(ID,start=4))

# Prepare a matrix of FC and p-values from each comparison in each cell line ------------
args <- c(number = nrow(mat), resort.by = "logFC")

temp <- data.frame("Genes" = rownames(topTable(fit_bayes, number=nrow(mat))))
for (i in seq.int(1,6,1)) {
  temp <- left_join(temp,(topTable(fit_bayes,coef=i,substitute(args))) |> dplyr::select(logFC) %>% rownames_to_column(var = "Genes"), by="Genes")
  temp <- left_join(temp,(topTable(fit_bayes,coef=i+6,substitute(args))) |> dplyr::select(logFC) %>% rownames_to_column(var = "Genes"), by="Genes")
  temp <- left_join(temp,(topTable(fit_bayes,coef=i+12,substitute(args))) |> dplyr::select(logFC) %>% rownames_to_column(var = "Genes"), by="Genes")
  colnames(temp)[(ncol(temp)-2):ncol(temp)] <- str_c(cell.line,curr.comp[i],"logFC", sep = ".")
}
assign("FC", temp)
FC <- column_to_rownames(FC, var="Genes")

# Prepare a matrix of means in each comparison in each cell line -------------------------------
mat <- as.data.frame(mat)
means.col <- colnames(mat) %>% str_sub(end=-3L) %>% unique()
gr.means <- data.frame(matrix(ncol=1,nrow=nrow(mat)))
for (i in means.col) {
 gr.means <- bind_cols(gr.means,mat %>% rowwise() %>% 
          transmute(mean(c_across(matches(paste0(i,"_[1234]"))), na.rm=TRUE)))
}
gr.means <- gr.means[-1]
colnames(gr.means) <- means.col
rownames(gr.means) <- rownames(mat)

# 5. Retrieve FC of the genes from each comparison in each cell line
term.genes.vec <- term.desc$Genes %>% flatten() %>% unlist() %>% unique()
FC_wp <- FC %>% filter(rownames(FC) %in% term.genes.vec)
#FC_kegg <- FC %>% filter(rownames(FC) %in% term.genes.vec)
#FC_bp <- FC %>% filter(rownames(FC) %in% term.genes.vec)
#FC_cc <- FC %>% filter(rownames(FC) %in% term.genes.vec)
#FC_mf <- FC %>% filter(rownames(FC) %in% term.genes.vec)

# 7. Draw heatmap
# HM loop for all terms --------------------------------------------------------------
terms.hm(FC.matrix = FC_wp, term.descr = term.desc, orderGroupWise = FALSE, save = TRUE, path = "figures/allruns/final_quant/term_overlap/WP/")
terms.hm(mean.matrix = gr.means, term.descr = term.desc, orderGroupWise = FALSE, save = TRUE, path = "figures/allruns/final_quant/term_overlap/WP/")

# KEGG panc cancer gene set tested ---------------------------------------------
# Parse GMT file for KEGG pancreatic cancer hallmark genes

#temp <- read.gmt(gmtfile = "data/KEGG_PANCREATIC_CANCER.v2024.1.Hs.gmt") # this one doesn't work for now
panc.canc <- clusterProfiler::read.gmt.wp(gmtfile = "data/GMTs/KEGG_PANCREATIC_CANCER.v2024.1.Hs.gmt")

files <- list.files("data/GMTs/", 
                    full.names = TRUE,
                    pattern = "VANGURP")

norm.panc <- lapply(files,clusterProfiler::read.gmt.wp)
norm.panc <- data.table::rbindlist(norm.panc)

pheatmap::pheatmap(gr.means %>% filter(rownames(gr.means) %in% panc.canc$gene),
                   cluster_cols = FALSE,
                   gaps_col = c(5,10,15),
                   cluster_rows = TRUE,
                   angle_col = 45,
                   border_color = "black",
                   col=colorRampPalette(c("cornflowerblue","white","red"))(100))
