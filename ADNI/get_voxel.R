library(tidyverse)
library(data.table)

sc <- function(dat, variable, variable.index, threshold){
  pvalue.adjust <- p.adjust(dat[, variable.index], method = 'fdr')
  index.dat <- which(pvalue.adjust<threshold)
  rows <- length(index.dat)
  if (rows==0){
    message <- paste0('if we choose threshold=', threshold, ' for variable ',
                      variable, ', then fdr adjusted p-value has no rows!')
    stop(message)
  }
  if (rows<100){
    message <- paste0('if we choose threshold=', threshold, ' for variable ',
                      variable, ', then fdr adjusted p-value has only ', rows, ' rows!')
    warning(message)
  } else{
    message <- paste0('if we choose threshold=', threshold, ' for variable ',
                      variable, ', then fdr adjusted p-value has ', rows, ' rows!\n')
    cat(message)
  }
  return(index.dat)
}

S <- 1843303
M <- 647
index <- NULL

surv.p <- data.frame(Value=rep(0, S), StdErr=rep(0, S), t_value=rep(0, S), p_value=rep(0, S))
  
## Don't run this for loop 
for (i in 1:M){
  start <- (i-1)*S/M+1; end <- i*S/M
  filename <- paste0("result/voxel_selection/", i, ".RData")
  load(filename)
  surv.p[start:end, ] <- p_mat
}
  
save(list = 'surv.p', file= 'RData/surv.p.RData')
#sc(surv.p, 'Survival', 4, threshold = 0.05)
#sc(surv.p, 'Survival', 4, threshold = 0.025)
#sc(surv.p, 'Survival', 4, threshold = 0.01)

index <- sc(surv.p, 'Survival', 4, threshold = 0.01)
save(list = 'index', file = 'RData/p_value_index.RData')
len.index = length(index)


mri.dir <- "/pine/scr/h/a/haotian/MMFPCA_MRI/ADNI_MRI/"
mri.files <- list.files(mri.dir, full.names = T)

surv = read.csv('dataset/ADNI_surv.csv')
N <- nrow(surv) ## number of subjects

Surv.mri.mat <-  matrix(NA, nrow = N, ncol = len.index)

for (i in 1:N){  
  tmp.dat <- fread(mri.files[i])
  Surv.mri.mat[i, ] <- tmp.dat$value[index]
}

write.csv(Surv.mri.mat, "dataset/Surv.mri.dat.csv", row.names = F)