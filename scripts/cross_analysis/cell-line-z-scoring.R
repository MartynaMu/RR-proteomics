# Normalization ---------------------------------------------------
#qnorm <- limma::normalizeQuantiles(as.matrix(df))
qnorm <- limma::normalizeMedianValues(as.matrix(df))

# Subtract means-----------------------------------------------------
# identify columns of each cell line
panc.ind <- colnames(qnorm) %>% str_which(pattern="PANC")
miapaca.ind <- colnames(qnorm) %>% str_which(pattern="MIAPACA")
cfpac.ind <- colnames(qnorm) %>% str_which(pattern="CFPAC")

# calculate means and sd of each cell line
temp <- qnorm %>% 
  as.data.frame() %>% 
  rowwise() %>% 
  mutate(PANC.mean = mean(c_across(all_of(panc.ind)),na.rm=TRUE),
         PANC.sd = sd(c_across(all_of(panc.ind)), na.rm = TRUE),
         MIAPACA.mean = mean(c_across(all_of(miapaca.ind)), na.rm = TRUE),
         MIAPACA.sd = sd(c_across(all_of(miapaca.ind)), na.rm = TRUE),
         CFPAC.mean = mean(c_across(all_of(cfpac.ind)), na.rm = TRUE),
         CFPAC.sd = sd(c_across(all_of(cfpac.ind)), na.rm = TRUE))

# subtract the means of cell lines from each
zscored <- mutate(temp, across(.cols = all_of(panc.ind), .fns = function(x) ((x - PANC.mean)/PANC.sd)),
              across(.cols = all_of(miapaca.ind), .fns = function(x) ((x - MIAPACA.mean)/MIAPACA.sd)),
              across(.cols = all_of(cfpac.ind), .fns = function(x) ((x - CFPAC.mean)/CFPAC.sd)))
zscored <- zscored %>% as.data.frame() %>% drop_na()
zscored <- zscored[1:48]
rownames(zscored) <- df %>% drop_na() %>% rownames()