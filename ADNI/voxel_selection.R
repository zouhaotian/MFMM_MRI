setwd('/pine/scr/h/a/haotian/MMFPCA_MRI/ADNI/')
library(tidyverse)
library(data.table)
library(survival)

# index = 1

mri.dir <- "/pine/scr/h/a/haotian/MMFPCA_MRI/ADNI_MRI/"
mri.files <- list.files(mri.dir, full.names = T)

## Read in the dataset ##
surv = read.csv('dataset/ADNI_surv.csv')
N <- nrow(surv) ## number of subjects
S <- 1843303 ## number of voxels 
M <- 647 ## M pieces, each with S/M voxels ##

mat <- matrix(rep(0, N*S/M), nrow=N)  ## records the volume data from start to end

start <- (index-1)*(S/M) + 1
end <- index*(S/M) 
for (i in 1:N){  
  tmp.mat <- fread(mri.files[i])
  mat[i, ] <- tmp.mat$value[start:end]
}

p_mat <- matrix(NA, nrow = S/M, ncol = 4)
tmp.voxels <- S/M

for (s in 1:tmp.voxels){
    tmp.volume <- mat[, s]
    cox.obj <- coxph(Surv(surv_time, status) ~ Age + Gender + Education + APOE4 + tmp.volume, 
                     method = 'breslow', data = surv)
    p_mat[s, ] <- summary(cox.obj)$coefficients[5, c(1, 3, 4, 5)]
} 

filename <- paste0('result/voxel_selection/', index, '.RData')
save(list = 'p_mat', file = filename)
