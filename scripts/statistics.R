library(tidyverse)

stats <- data.frame(x2D_mean = rowMeans(temp[,1:4]), 
                 xPDX_2D_mean = rowMeans(temp[,13:15]))
stats <- stats %>% mutate(fc_2D_PDX2D = x2D_mean-xPDX_2D_mean)

temp <- rownames_to_column(as.data.frame(mat), "Protein.name")
temp <- as.data.frame(temp)
temp <- pivot_longer(temp, 2:19, names_to = "Sample", values_to = "Intensity")
temp2 <- temp |> mutate(
  Group1 = case_when(
    Sample %in% c("X2D1", "X2D2", "X2D3", "X2D4") ~ "2D",
    Sample %in% c("X3Y1", "X3Y2", "X3Y3", "X3Y4", "X3O1", "X3O2", "X3O3", "X3O4") ~ "3D",
    Sample %in% c("X2A","X2B","X2C","X3A","X3B","X3C") ~ "PDX"),
  Group2 = case_when(
    Sample %in% c("X2D1", "X2D2", "X2D3", "X2D4") ~ "2D",
    Sample %in% c("X3O1", "X3O2", "X3O3", "X3O4") ~ "3D old",
    Sample %in% c("X3Y1", "X3Y2", "X3Y3", "X3Y4") ~ "3D young",
    Sample %in% c("X2A","X2B","X2C") ~ "PDX 2D",
    Sample %in% c("X3A","X3B","X3C") ~ "PDX 3D"
  )
)
x <- t.test(Intensity ~ Group1, data = temp2[temp2$Group1 %in% c("2D", "3D"),])
temp3 <- rowwise(mutate(temp2, pval = x$p.value))

