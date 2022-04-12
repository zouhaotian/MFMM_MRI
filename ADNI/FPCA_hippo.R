## Hippocampus FPCA ##
library(data.table)
library(tidyverse)
library(splines)
library(mgcv)
library(MBA)
library(RSpectra)

source('source_code/f0.R')
MRI.dat <- fread("dataset/hippo.mri.dat.csv") %>% as.matrix()

FPCA2 <- function(dat){
  S <- ncol(dat)
  obsgrid <- 0:(S-1)
  by <- 1
  FPCA.obj <- FPCA_hd(dat, obsgrid, pve = 0.95, by = by, L = 8)
  return(FPCA.obj)
}
FPCA.res <- FPCA2(MRI.dat)

save(list = 'FPCA.res', 
     file = 'RData/FPCA_Hippo_result.RData')