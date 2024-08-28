# Use to plot and save heatmaps as png for overlapped enriched terms
terms.hm <- function(FC.matrix = NULL, mean.matrix = NULL, term.descr, orderGroupWise = FALSE, save = TRUE, path = NULL) {
 
  terms <- term.descr$ID
  
  if (!is.null(FC.matrix)) {
    matrix <- FC.matrix
    plot.desc <- "Log2FC in each comparison"
    file.desc <- ".hm.fc.png"
    path2 <- "fc/"
    gaps_col <- c(3,6,9,12,15)
  } else if (!is.null(mean.matrix)) {
    matrix <- mean.matrix
    plot.desc <- "Mean group intensity"
    file.desc <- ".hm.mean.png"
    path2 <- "means/"
    gaps_col <- c(5,10,15)
    if (orderGroupWise==TRUE) {
      # order colnames group-wise
      order <- colnames(mean.matrix) %>% str_replace_all(pattern = "MIAPACA|PANC|CFPAC", replacement = "") %>% unique()
      order[1] <- str_c(order[1],"$")
      
      temp <- c()
      for (i in order) {
        temp <- c(temp, str_which(colnames(mean.matrix), pattern = i))
      }
      mean.matrix <- relocate(mean.matrix, all_of(temp))
      gaps_col <- c(3,6,9,12,15)
    }
  }
  
  for (i in terms) {
    term.desc <- term.descr %>% filter(ID == i) %>% pull(Description)
    term.genes <- term.descr %>% filter(ID == i) %>% pull(Genes)
    names(term.genes) <- i
    matrix.to.use <- matrix %>% filter(rownames(matrix) %in% term.genes[[i]])
    if (length(term.genes[[i]]) > 20) {
      p <- pheatmap::pheatmap(matrix.to.use,
                              cluster_cols = FALSE,
                              gaps_col = c(5,10,15),
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
                              gaps_col = c(5,10,15),
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
