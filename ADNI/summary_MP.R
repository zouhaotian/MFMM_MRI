library(tidyverse)
library(xtable)

models.name = c('posterior_M1_whole', 'posterior_M1', 'prior_M1', 
                'posterior_M2_whole', 'posterior_M2', 'prior_M2',  
                'posterior_M3_2S', 'posterior_M3_MJM_MRI', 
                'posterior_M3_bCox', 'posterior_M3_non')
models.label = c('Model 1 - whole', 'Model 1 - hippo', 'Model 1', 
                 'Model 2 - whole', 'Model 2 - hippo', 'Model 2', 
                 'Model 3 - 2S', 'Model 3 - MJM-MRI', 
                 'Model 3 - bCox', 'Model 3 - NM')

summary <- matrix(NA, nrow = length(models.name), ncol = 5)
for (m in 1:length(models.name)){
  fname = paste0('RData2/summary_', models.name[m], '_MP.RData')
  load(fname)
  summary[m, 1] <- l$EAIC
  summary[m, 2] <- l$EBIC
  summary[m, 3] <- l$looic
  summary[m, 4] <- l$waic
}

summary[, 5] <- c(11941.4, 9409.75, 9172.05, 
                  21159.1, 28946.6, 24165.1, 
                  17833.8, 41089.1, 20740.2, 3515.93)/3600

rownames(summary) <- models.label
colnames(summary) <- c('EAIC', 'EBIC', 'LOOIC', 'WAIC', 'Time (hrs)')

write.csv(summary, file = 'result/summary_MP.csv', row.names = T)

summ = read.csv('result/summary_MP.csv')
print(xtable(summ, digits = c(0, rep(0, 5), 1)), include.rownames = F)

