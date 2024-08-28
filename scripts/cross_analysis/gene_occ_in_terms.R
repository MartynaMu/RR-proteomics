cell.line <- c("miapaca", "panc1", "cfpac")

gene_occ <- data.frame(matrix(ncol=1))
for (i in wp_list) {
  # count terms that every gene was taken account into
  curr_gene_occ <- data.frame(matrix(ncol=2))
  colnames(curr_gene_occ) <- c("Gene", "TermCount")
  
  for (g in rownames(mat)) {
    temp <- i@result$core_enrichment |> str_count(pattern = g) %>% sum()
    curr_gene_occ <- add_row(curr_gene_occ, Gene = g, TermCount = temp)
  }
  curr_gene_occ <- curr_gene_occ %>% drop_na() %>% arrange(desc(TermCount))
  print(head(curr_gene_occ))
  
 # assign(paste("gene_occ",names(coefs[i]),sep="."), gene_occ)
}


test <- full_join(gene_occ.cfpac.2dv3dy, gene_occ.miapaca.2dv3dy, by="Gene") |> full_join(gene_occ.panc1.2dv3dy, by="Gene")

test1 <- test %>% mutate_at(vars(TermCount.x:TermCount), rank, na.last = "keep")

test1 |> pivot_longer(2:4, names_to = "CellLine", values_to = "TermCount") %>%
ggplot(aes(x=Gene,y=TermCount, color = CellLine, group=CellLine))+
geom_line()

