# Evaluate min and max of term lengths
temp <- data.frame(matrix(ncol=1))
colnames(temp) <- "Length"
all.lists <- list(bp_list, cc_list, mf_list, wp_list)

for (i in all.lists) {
  for (e in 1:18) {
    curr <- i[[e]]@result 
    curr <- curr %>% rowwise() %>% mutate(core_enrichment= str_split(core_enrichment, pattern = "/") |> map(unique))
    curr <- curr %>% rowwise() %>% transmute(Length=length(core_enrichment))
    temp <- rbind(temp,curr)
  }
}
temp <- drop_na(temp)
temp %>% ggplot(aes(x=Length))+ geom_histogram()
min(temp$Length)
max(temp$Length)
temp1 <- data.frame(A=1:83, B=sample.int(83,size = 83))
pheatmap::pheatmap(temp1[1:80,], cluster_cols = FALSE, fontsize_row = 7)



# Function
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
    term.genes <- term.descr %>% filter(ID == i) %>% pull(Genes) %>% unlist()
    # determine figure height
    if(length(term.genes) < 20){
      height <- 5
    } else if(length(term.genes) > 20 && length(term.genes) <= 40) {
      height <- 6
    } else if(length(term.genes) > 40 && length(term.genes) <= 60) {
      height <- 7
    } else if(length(term.genes) > 60 && length(term.genes) <= 90) {
      height <- 8
    }
    
    # Means grouped cell.line-wise
    cat("  \n##", sprintf("%s (%s)",curr.term.desc, i), "  \n")
    cat("  \n## {.tabset}  \n")
    cat("  \n### Means cell line-wise  \n")
    if (length(term.genes) > 20) {
      p <- mean.matrix %>% filter(rownames(mean.matrix) %in% term.genes) |>
        pheatmap::pheatmap(cluster_cols = FALSE,
                           cluster_rows = TRUE,
                           angle_col = 45,
                           fontsize_row = 7,
                           border_color = "black",
                           gaps_col = c(5,10,15),
                           col=colorRampPalette(c("cornflowerblue","white","red"))(100))
      cat("  \n")
    } else if (length(term.genes) <= 20) {
      p <- mean.matrix %>% filter(rownames(mean.matrix) %in% term.genes) |>
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
    }
    # Means grouped experimental group-wise
    cat("  \n### Means exp.group-wise  \n")
    if (length(term.genes) > 20) {
      p <- mean.matrix.gr %>% filter(rownames(mean.matrix.gr) %in% term.genes) |>
        pheatmap::pheatmap(cluster_cols = FALSE,
                           cluster_rows = TRUE,
                           angle_col = 45,
                           fontsize_row = 7,
                           border_color = "black",
                           gaps_col = c(3,6,9,12,15),
                           col=colorRampPalette(c("cornflowerblue","white","red"))(100))
      cat("  \n")
    } else if (length(term.genes) <= 20) {
      p <- mean.matrix.gr %>% filter(rownames(mean.matrix.gr) %in% term.genes) |>
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
    }
    # Log2FC matrix
    cat("  \n### Log2FC  \n")
    if (length(term.genes) > 20) {
      p <- FC.matrix %>% filter(rownames(FC.matrix) %in% term.genes) |>
        pheatmap::pheatmap(cluster_cols = FALSE,
                           cluster_rows = TRUE,
                           angle_col = 45,
                           fontsize_row = 7,
                           border_color = "black",
                           gaps_col = c(3,6,9,12,15),
                           col=colorRampPalette(c("cornflowerblue","white","red"))(100))
      cat("  \n")
    } else if (length(term.genes) <= 20) {
      p <- FC.matrix %>% filter(rownames(FC.matrix) %in% term.genes) |>
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
    }
  }
}
