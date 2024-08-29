cell.line <- c("miapaca", "panc1", "cfpac")

gene_occ <- list()
for (i in seq.int(1,18,1)) {
  # count terms that every gene was taken account into
  curr_gene_occ <- data.frame(matrix(ncol=3))
  colnames(curr_gene_occ) <- c("Gene", "TermCount", "Direction")
  
  for (g in rownames(mat)) {
    temp.up <- wp_list[[i]]@result %>% filter(NES > 0) %>% pull(core_enrichment) |> str_count(pattern = g) %>% sum()
    curr_gene_occ <- add_row(curr_gene_occ, Gene = g, TermCount = temp.up, Direction = "Up")
    temp.down <- wp_list[[i]]@result %>% filter(NES < 0) %>% pull(core_enrichment) |> str_count(pattern = g) %>% sum()
    curr_gene_occ <- add_row(curr_gene_occ, Gene = g, TermCount = temp.down, Direction = "Down")
  }
  curr_gene_occ <- curr_gene_occ %>% drop_na() %>% arrange(desc(TermCount))
  print(head(curr_gene_occ))
  gene_occ[[i]] <- curr_gene_occ
}

names(gene_occ) <- names(coefs)
gene_occ <- enframe(gene_occ) %>% unnest()
gene_occ[gene_occ==0]<- NA
gene_occ <- drop_na(gene_occ)
gene_occ <- separate(gene_occ, col=name,into=c("Line", "Comp"), sep="\\.")
comp_gene_occ <- gene_occ %>% 
  group_by(Comp, Gene, Direction) %>% 
  summarise(CompCount = sum(TermCount)) %>% 
  arrange(desc(CompCount))
sum_gene_occ <- comp_gene_occ %>% 
  group_by(Gene, Direction) %>% 
  summarise(Count = sum(CompCount)) %>% 
  arrange(desc(Count))

sum_gene_occ %>% 
  filter(Count>30) |> 
  ggbarplot(x = "Gene", y="Count", 
            fill = "Direction",
            palette = c("steelblue","firebrick"),
            label=TRUE,
            xlab="Gene symbol",
            ylab="Gene count",
            title = "Top gene count in WikiPathways GSEA-enriched terms in all comparisons") |>
  facet(facet.by = "Direction", nrow=2)+
  font("x.text", size = 8)+
  scale_y_continuous(limits = c(0,120))

filtr2 <- sum_gene_occ %>% filter(Count > 30) %>% pull(Gene)
comp_gene_occ %>% 
  mutate(CompCount = case_when(Direction == "Down" ~ CompCount*-(1), .default = CompCount)) |>
  filter(Gene %in% filtr2) |>
  ggbarplot(x="Gene", y="CompCount",
            fill = "Direction",
            palette = c("steelblue","firebrick"),
            rotate=TRUE,
            xlab = "Gene symbol",
            ylab = "Gene count",
            title = "Top gene count in WikiPathways GSEA-enriched terms in each comparison") |>
  facet(facet.by = c("Comp"), ncol = 2,panel.labs.font = list(size=12), 
        scales="free_y",
        panel.labs.background = list(color="white", fill="white"))+
  scale_y_continuous(limits=c(-45,45))+
  geom_text(aes(label=CompCount), hjust = "outward", size=3.5)+
  font("xy.text", size = 8)


