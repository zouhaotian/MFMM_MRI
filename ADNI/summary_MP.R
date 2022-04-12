library(tidyverse)
library(xtable)

models.name = c('prior_M1', 'posterior_M1', 'posterior_M1_whole', 
                'prior_M2', 'posterior_M2', 'posterior_M2_whole')
models.label = c('M1_MFMM', 'M1_hippo', 'M1_whole', 
                 'M2_MFMM', 'M2_hippo', 'M2_whole')

summary <- matrix(NA, nrow = length(models.name), ncol = 4)
for (m in 1:length(models.name)){
  fname = paste0('RData/summary/summary_', models.name[m], '_MP.RData')
  load(fname)
  summary[m, 1] <- l$EAIC
  summary[m, 2] <- l$EBIC
  summary[m, 3] <- l$looic
  summary[m, 4] <- l$waic
}

rownames(summary) <- models.label
colnames(summary) <- c('EAIC', 'EBIC', 'LOOIC', 'WAIC')

write.csv(summary, file = 'result/summary_MP.csv')

summ = read.csv('result/summary_MP.csv')
xtable(summ[, -1], digits = c(0, rep(3, 4)))
