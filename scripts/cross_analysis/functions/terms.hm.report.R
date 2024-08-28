terms.hm.report <- function(FC.matrix, mean.matrix, term.descr) {
  
  terms <- term.descr$ID
  
  # order colnames group-wise
  order <- colnames(mean.matrix) %>% str_replace_all(pattern = "MIAPACA|PANC|CFPAC", replacement = "") %>% unique()
  order[1] <- str_c(order[1],"$")
  
  temp <- c()
  for (i in order) {
    temp <- c(temp, str_which(colnames(mean.matrix), pattern = i))
  }
  mean.matrix.gr <- relocate(mean.matrix, all_of(temp))
  
  for (i in terms) {
    curr.term.desc <- term.descr %>% filter(ID == i) %>% pull(Description)
    term.genes <- term.descr %>% filter(ID == i) %>% pull(Genes)
    names(term.genes) <- i
    # Means grouped cell.line-wise
    cat("  \n##", sprintf("%s (%s)",curr.term.desc, i), "  \n")
    cat("  \n## {.tabset}  \n")
    cat("  \n### Means cell line-wise  \n")
    if (length(term.genes[[i]]) > 20) {
      p <- mean.matrix %>% filter(rownames(mean.matrix) %in% term.genes[[i]]) |>
        pheatmap::pheatmap(cluster_cols = FALSE,
                           cluster_rows = TRUE,
                           angle_col = 45,
                           fontsize_row = 7,
                           border_color = "black",
                           gaps_col = c(5,10,15),
                           col=colorRampPalette(c("cornflowerblue","white","red"))(100))
      cat("  \n")
      #p
      cat("  \n")
    } else if (length(term.genes[[i]]) <= 20) {
      p <- mean.matrix %>% filter(rownames(mean.matrix) %in% term.genes[[i]]) |>
        pheatmap::pheatmap(cluster_cols = FALSE,
                           cluster_rows = TRUE,
                           angle_col = 45,
                           display_numbers = TRUE,
                           fontsize_number = 10,
                           fontsize_row = 7,
                           number_color = "black",
                           border_color = "black",
                           gaps_col = c(5,10,15),
                           col=colorRampPalette(c("cornflowerblue","white","red"))(100))
      cat("  \n")
      #p
      cat("  \n")
    }
    # Means grouped experimental group-wise
    cat("  \n### Means exp.group-wise  \n")
    if (length(term.genes[[i]]) > 20) {
      p <- mean.matrix.gr %>% filter(rownames(mean.matrix.gr) %in% term.genes[[i]]) |>
        pheatmap::pheatmap(cluster_cols = FALSE,
                           cluster_rows = TRUE,
                           angle_col = 45,
                           fontsize_row = 7,
                           border_color = "black",
                           gaps_col = c(3,6,9,12,15),
                           col=colorRampPalette(c("cornflowerblue","white","red"))(100))
      cat("  \n")
      #p
      cat("  \n")
    } else if (length(term.genes[[i]]) <= 20) {
      p <- mean.matrix.gr %>% filter(rownames(mean.matrix.gr) %in% term.genes[[i]]) |>
        pheatmap::pheatmap(cluster_cols = FALSE,
                           cluster_rows = TRUE,
                           angle_col = 45,
                           display_numbers = TRUE,
                           fontsize_number = 10,
                           fontsize_row = 7,
                           number_color = "black",
                           border_color = "black",
                           gaps_col = c(3,6,9,12,15),
                           col=colorRampPalette(c("cornflowerblue","white","red"))(100))
      cat("  \n")
      #p
      cat("  \n")
    }
    # Log2FC matrix
    cat("  \n### Log2FC  \n")
    if (length(term.genes[[i]]) > 20) {
      p <- FC.matrix %>% filter(rownames(FC.matrix) %in% term.genes[[i]]) |>
        pheatmap::pheatmap(cluster_cols = FALSE,
                           cluster_rows = TRUE,
                           angle_col = 45,
                           fontsize_row = 7,
                           border_color = "black",
                           gaps_col = c(3,6,9,12,15),
                           col=colorRampPalette(c("cornflowerblue","white","red"))(100))
      cat("  \n")
      #p
      cat("  \n")
    } else if (length(term.genes[[i]]) <= 20) {
      p <- FC.matrix %>% filter(rownames(FC.matrix) %in% term.genes[[i]]) |>
        pheatmap::pheatmap(cluster_cols = FALSE,
                           cluster_rows = TRUE,
                           angle_col = 45,
                           display_numbers = TRUE,
                           fontsize_number = 10,
                           fontsize_row = 7,
                           number_color = "black",
                           border_color = "black",
                           gaps_col = c(3,6,9,12,15),
                           col=colorRampPalette(c("cornflowerblue","white","red"))(100))
      cat("  \n")
      #p
      cat("  \n")
    }
  }
}

