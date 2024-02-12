

files <- list.files(path="data/gsego_results/gsego_BP_results/unfiltered/notsimplified", pattern=".tsv", full.names = TRUE)
one_file <- lapply(files,read_tsv)

temp <- do.call(rbind, one_file)


term = "negative regulation of cell cycle process"
duplicates <- unique(temp$Description[duplicated(temp$Description)])
str_flatten(temp$core_enrichment[temp$Description == ])
filtr <- str_split_1(temp$core_enrichment[temp$Description == term], pattern = "/")

pheatmap::pheatmap(qnorm[qnorm$Gene %in% filtr,1:18],
                   scale = "row",
                   show_rownames = TRUE,
                   clustering_distance_cols = "euclidean",
                   clustering_distance_rows = "correlation",
                   annotation_col = annot_col,
                   annotation_colors = annot_colors,
                   main = term)