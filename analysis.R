library(ggplot2)
library(cowplot)
library(dplyr)
library(kableExtra)
list.files(pattern = '.rds$')
date.file <- '2025_02_03'
time.file <- '1459'

# all-in-one analysis
fileName <- paste0(paste('simulationRes', date.file, time.file, 'summary', sep = '_'),'.rds')
res <- readRDS(fileName)

# accuracy table
accuaracy.tb <- lapply(res, function(res.) res.$accuracy.tb)
accuaracy.tb <- do.call(rbind, accuaracy.tb)
accuaracy.tb.sparse <- subset(accuaracy.tb, sparse.lv == 'sparse')
accuaracy.tb.neither <- subset(accuaracy.tb, sparse.lv == 'neither')
accuaracy.tb.dense <- subset(accuaracy.tb, sparse.lv == 'dense')
# scenarios - copied from simulation file.
eigenVals <- list(lam.1 = c(9, 4, 1), lam.2 = c(36, 25, 16, 9, 4, 1))
maxK <- c(6, 12) # maximum number of fPCs under consideration
sig2ma <- c(0.1, 0.5, 1, 4) # variance for random errors
n <- list(n.1 = c(100, 200), n.2 = c(200, 300)) # sample size
numObs <- c(11, 26, 51) # number of observations per curve
M <- 500 # simulation size
sparse.lv <- c('sparse', 'neither', 'dense')
scenarios <- lapply(1:length(eigenVals), function(i.) {
  expand.grid(n = n[[i.]], sparse.lv = sparse.lv, eigenvals = eigenVals[i.],
              sig2ma = sig2ma, isGau = c(TRUE, FALSE), isRegular = c(TRUE, FALSE),
              stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
})

scenarios <- do.call(rbind, scenarios)
# if the data is sparse/neither, we do not consider regular collected data
scenarios <- subset(scenarios, !(sparse.lv %in% c('sparse', 'neither') & isRegular == TRUE))
#scenarios$isRegular <- ifelse(scenarios$sparse.lv == 'dense', TRUE, FALSE)
scenarios$maxK <- ifelse(names(scenarios$eigenvals) == 'lam.1', maxK[1], maxK[2])
scenarios <- cbind(scenario = 1:nrow(scenarios), scenarios)
row.names(scenarios) <- NULL

# contingency table - for each scenario only
# case 81, 82, 83, 84
case <- 83
contingency.tb <- res[[case]]$occurence.tb
(setting <- scenarios[case,])
contingency.tb

# output table into latex
selected.cols <- c("scenario", "n", "sparse.lv", "eigenvals", "sig2ma", "isGau",
                   "isRegular", "k.prop.wn1", "k.prop.wn2", "k.prop.wn3", "k.prop.wn4",
                   "k.prop.won", "FY.AIC", "FY.BIC", "LYH.BIC", "LYH.AIC")
accuaracy.tb <- accuaracy.tb[, selected.cols]


#' table generating 
library(dplyr)
library(tidyr)
library(tidyverse)
methods_names <- c("k.prop.wn1", "k.prop.wn2", "k.prop.wn3", "k.prop.wn4",
                   "k.prop.won", "FY.AIC", "FY.BIC", "LYH.BIC", "LYH.AIC")
table_generator <- function(raw.tb, eigenSetting, is.Gau, is.Regular, sparseLv, method_names){
  if(sparseLv == 'dense'){
    temp.tb <- subset(raw.tb, names(eigenvals) == eigenSetting & isGau == is.Gau
                      & sparse.lv == sparseLv & isRegular == is.Regular)
  } else{
    temp.tb <- subset(raw.tb, names(eigenvals) == eigenSetting & isGau == is.Gau
                      & sparse.lv == sparseLv)
  }
  
  temp.ltb <- pivot_longer(data = temp.tb, names_to = 'method',
                           values_to = 'accuracy', cols = 8:16,
                           cols_vary = "fastest")
  temp.ltb <- temp.ltb[, -which(names(temp.ltb) == 'scenario')]
  temp.wtb <- 
    pivot_wider(data = temp.ltb, names_from = sig2ma,
                names_prefix = "sig2ma.", values_from = accuracy)
  n.seq <- sort(unique(temp.wtb$n))
  mat.l <- lapply(n.seq, function(n.){
    temp.df <- subset(temp.wtb, n == n.,
                      select = which(str_detect(names(temp.wtb), '^sig2m')))
    temp.mat <- as.matrix(temp.df)
    #colnames(temp.mat) <- paste0(colnames(temp.mat), 'n.',n.)
    row.names(temp.mat) <- method_names
    return(temp.mat)
  })
  final.tb <- data.frame(do.call(cbind, mat.l))
  final.tb <- cbind(final.tb, unique(temp.wtb[, 2:5]))
  return(final.tb)
}
###------------------------ Gaussian -------------------------------
#' table 1 - simple + Gau + sparse/neither
accuaracy.simple.GauAndSparse.wtb <- 
  table_generator(raw.tb = accuaracy.tb, eigenSetting = 'lam.1', is.Gau = TRUE,
                  is.Regular = FALSE, sparseLv = 'sparse', method_names = methods_names)
accuaracy.simple.GauAndSparse.wtb %>% kable(, format = 'latex')

accuaracy.simple.GauAndNeither.wtb <- 
  table_generator(raw.tb = accuaracy.tb, eigenSetting = 'lam.1', is.Gau = TRUE,
                  is.Regular = FALSE, sparseLv = 'neither', method_names = methods_names)
accuaracy.simple.GauAndNeither.wtb %>% kable(, format = 'latex')

#' table 2 - simple + Gau + dense + regular/irregular
accuaracy.simple.GauAndDenseAndIrr.wtb <- 
  table_generator(raw.tb = accuaracy.tb, eigenSetting = 'lam.1', is.Gau = TRUE,
                  is.Regular = FALSE, sparseLv = 'dense', method_names = methods_names)
accuaracy.simple.GauAndDenseAndIrr.wtb %>% kable(, format = 'latex')

accuaracy.simple.GauAndDenseAndReg.wtb <- 
  table_generator(raw.tb = accuaracy.tb, eigenSetting = 'lam.1', is.Gau = TRUE,
                  is.Regular = TRUE, sparseLv = 'dense', method_names = methods_names)
accuaracy.simple.GauAndDenseAndReg.wtb %>% kable(, format = 'latex', digits = 3)

#' table 3 - Complex + Gau + sparse/neither
accuaracy.complex.GauAndSparse.wtb <- 
  table_generator(raw.tb = accuaracy.tb, eigenSetting = 'lam.2', is.Gau = TRUE,
                  is.Regular = FALSE, sparseLv = 'sparse', method_names = methods_names)
accuaracy.complex.GauAndSparse.wtb %>% kable(, format = 'latex')

accuaracy.complex.GauAndNeither.wtb <- 
  table_generator(raw.tb = accuaracy.tb, eigenSetting = 'lam.2', is.Gau = TRUE,
                  is.Regular = FALSE, sparseLv = 'neither', method_names = methods_names)
accuaracy.complex.GauAndNeither.wtb %>% kable(, format = 'latex')

#' table 4 - Complex + Gau + dense + regular/irregular
accuaracy.complex.GauAndDenseAndIrr.wtb <- 
  table_generator(raw.tb = accuaracy.tb, eigenSetting = 'lam.2', is.Gau = TRUE,
                  is.Regular = FALSE, sparseLv = 'dense', method_names = methods_names)
accuaracy.complex.GauAndDenseAndIrr.wtb %>% kable(, format = 'latex')

accuaracy.complex.GauAndDenseAndReg.wtb <- 
  table_generator(raw.tb = accuaracy.tb, eigenSetting = 'lam.2', is.Gau = TRUE,
                  is.Regular = TRUE, sparseLv = 'dense', method_names = methods_names)
accuaracy.complex.GauAndDenseAndReg.wtb %>% kable(, format = 'latex', digits = 3)
###------------------------ non-Gaussian -------------------------------
#' table 5 - Simple + nonGau + sparse/neither
accuaracy.simple.nonGauAndSparse.wtb <- 
  table_generator(raw.tb = accuaracy.tb, eigenSetting = 'lam.1', is.Gau = FALSE,
                  is.Regular = FALSE, sparseLv = 'sparse', method_names = methods_names)
accuaracy.simple.nonGauAndSparse.wtb %>% kable(, format = 'latex')

accuaracy.simple.nonGauAndNeither.wtb <- 
  table_generator(raw.tb = accuaracy.tb, eigenSetting = 'lam.1', is.Gau = FALSE,
                  is.Regular = FALSE, sparseLv = 'neither', method_names = methods_names)
accuaracy.simple.nonGauAndNeither.wtb %>% kable(, format = 'latex')

#' table 6 - Simple + nonGau + dense + reg/irreg
accuaracy.simple.nonGauAndDenseAndIrr.wtb <- 
  table_generator(raw.tb = accuaracy.tb, eigenSetting = 'lam.1', is.Gau = FALSE,
                  is.Regular = FALSE, sparseLv = 'dense', method_names = methods_names)
accuaracy.simple.nonGauAndDenseAndIrr.wtb %>% kable(, format = 'latex')

accuaracy.simple.nonGauAndDenseAndReg.wtb <- 
  table_generator(raw.tb = accuaracy.tb, eigenSetting = 'lam.1', is.Gau = FALSE,
                  is.Regular = TRUE, sparseLv = 'dense', method_names = methods_names)
accuaracy.simple.nonGauAndDenseAndReg.wtb %>% kable(, format = 'latex', digits = 3)

#' table 7 - Complex + nonGau + sparse/neither
accuaracy.complex.nonGauAndSparse.wtb <- 
  table_generator(raw.tb = accuaracy.tb, eigenSetting = 'lam.2', is.Gau = FALSE,
                  is.Regular = FALSE, sparseLv = 'sparse', method_names = methods_names)
accuaracy.complex.nonGauAndSparse.wtb %>% kable(, format = 'latex')

accuaracy.complex.nonGauAndNeither.wtb <- 
  table_generator(raw.tb = accuaracy.tb, eigenSetting = 'lam.2', is.Gau = FALSE,
                  is.Regular = FALSE, sparseLv = 'neither', method_names = methods_names)
accuaracy.complex.nonGauAndNeither.wtb %>% kable(, format = 'latex')

#' table 8 - Complex + nonGau + dense + reg/irreg
accuaracy.complex.nonGauAndDenseAndIrr.wtb <- 
  table_generator(raw.tb = accuaracy.tb, eigenSetting = 'lam.2', is.Gau = FALSE,
                  is.Regular = FALSE, sparseLv = 'dense', method_names = methods_names)
accuaracy.complex.nonGauAndDenseAndIrr.wtb %>% kable(, format = 'latex')

accuaracy.complex.nonGauAndDenseAndReg.wtb <- 
  table_generator(raw.tb = accuaracy.tb, eigenSetting = 'lam.2', is.Gau = FALSE,
                  is.Regular = TRUE, sparseLv = 'dense', method_names = methods_names)
accuaracy.complex.nonGauAndDenseAndReg.wtb %>% kable(, format = 'latex', digits = 3)
