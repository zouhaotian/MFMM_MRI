library(tidyverse)
library(data.table)
library(EveTemplate)

S <- 1843303
N <- 742
mri.dir <- "E:/MMFPCA_MRI/ADNI_MRI"
mri.file <- list.files(mri.dir, full.names = T)

eve_brain_mask <- readEve(what = 'Brain_Mask')
eve_labels <- readEveMap(type = 'II')
label_df <- getEveMapLabels("II")
t1 <- readEve("Brain")

df <- data.frame(integer_label = eve_labels[eve_brain_mask==1], index = c(1:S)) %>% 
  left_join(label_df, by = 'integer_label')
hippo.index <- df[which(df$structure=='hippocampus'), ]$index

hippo.mri.mat <- matrix(rep(0, N*length(hippo.index)), nrow=N)
for (i in 1:N){
  tmp.dat <- fread(mri.file[i])
  hippo.mri.mat[i, ] <- tmp.dat$value[hippo.index]
}

write.csv(hippo.mri.mat, "dataset/hippo.mri.dat.csv", row.names = F)