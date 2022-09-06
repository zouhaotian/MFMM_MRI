library(tidyverse)
library(xtable)

models.name = c('posterior_M2_whole_MP_P7', 'posterior_M2_whole_MP_P8', 
                'posterior_M2_whole_MP_P9', 'posterior_M2_whole_MP_P10', 
                'posterior_M2_whole_MP_P11')
models.label = c("M2_whole: P=7", "M2_whole: P=8", 
                 "M2_whole: P=9", "M2_whole: P=10", 
                 "M2_whole: P=11")

summary <- matrix(NA, nrow = length(models.name), ncol = 4)
for (m in 1:length(models.name)){
  fname = paste0('RData2/summary_', models.name[m], '.RData')
  load(fname)
  summary[m, 1] <- l$EAIC
  summary[m, 2] <- l$EBIC
  summary[m, 3] <- l$looic
  summary[m, 4] <- l$waic
}

rownames(summary) <- models.label
colnames(summary) <- c('EAIC', 'EBIC', 'LOOIC', 'WAIC')

write.csv(summary, file = 'result/summary_MP_sensitivity.csv')

summ = read.csv('result/summary_MP_sensitivity.csv')
xtable(summ[, -1], digits = c(0, rep(0, 4)))
