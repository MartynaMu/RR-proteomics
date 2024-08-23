terms.hm.report <- function(FC.matrix, mean.matrix, term.descr) {

  terms <- term.descr$ID
  for (i in terms) {
    curr.term.desc <- term.descr %>% filter(ID == i) %>% pull(Description)
    term.genes <- term.descr %>% filter(ID == i) %>% pull(Genes)
    names(term.genes) <- i
    
    cat("  \n## {.tabset}  \n")
    cat("  \n### Log2FC  \n")
    cat("  \n", sprintf("%s (%s)",curr.term.desc, i), "  \n")
    if (length(term.genes[[i]]) > 20) {
      p <- FC.matrix %>% filter(rownames(FC.matrix) %in% term.genes[[i]]) |>
        pheatmap::pheatmap(cluster_cols = FALSE,
                              gaps_col = c(3,6,9,12,15),
                              cluster_rows = TRUE,
                              angle_col = 45,
                              border_color = "black",
                              main = sprintf("Log2FC in each comparison - %s (%s)", curr.term.desc, i),
                              col=colorRampPalette(c("cornflowerblue","white","red"))(100))
      cat("  \n")
      #p
      cat("  \n")
    } else if (length(term.genes[[i]]) <= 20) {
      p <- FC.matrix %>% filter(rownames(FC.matrix) %in% term.genes[[i]]) |>
        pheatmap::pheatmap(cluster_cols = FALSE,
                              gaps_col = c(3,6,9,12,15),
                              cluster_rows = TRUE,
                              angle_col = 45,
                              display_numbers = TRUE,
                              fontsize_number = 10,
                              number_color = "black",
                              border_color = "black",
                              main = sprintf("Log2FC in each comparison - %s (%s)", curr.term.desc, i),
                              col=colorRampPalette(c("cornflowerblue","white","red"))(100))
      cat("  \n")
      #p
      cat("  \n")
    }
    cat("  \n### Means  \n")
    cat("  \n", sprintf("%s (%s)",curr.term.desc, i), "  \n")
    if (length(term.genes[[i]]) > 20) {
      p <- mean.matrix %>% filter(rownames(mean.matrix) %in% term.genes[[i]]) |>
        pheatmap::pheatmap(cluster_cols = FALSE,
                           gaps_col = c(3,6,9,12,15),
                           cluster_rows = TRUE,
                           angle_col = 45,
                           border_color = "black",
                           main = sprintf("Mean group intensity - %s (%s)", curr.term.desc, i),
                           col=colorRampPalette(c("cornflowerblue","white","red"))(100))
      cat("  \n")
      #p
      cat("  \n")
    } else if (length(term.genes[[i]]) <= 20) {
      p <- mean.matrix %>% filter(rownames(mean.matrix) %in% term.genes[[i]]) |>
        pheatmap::pheatmap(cluster_cols = FALSE,
                           gaps_col = c(3,6,9,12,15),
                           cluster_rows = TRUE,
                           angle_col = 45,
                           display_numbers = TRUE,
                           fontsize_number = 10,
                           number_color = "black",
                           border_color = "black",
                           main = sprintf("Mean group intensity - %s (%s)", curr.term.desc, i),
                           col=colorRampPalette(c("cornflowerblue","white","red"))(100))
      cat("  \n")
      #p
      cat("  \n")
    }
  }
}