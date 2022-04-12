library(tidyverse)


datadir <- 'RData/M2_whole_bs/'

result.beta_m <- NULL
index.res <- 0

for (i in 1:120){
  fname <- paste0(datadir, i, '.RData')
  if (file.exists(fname)){
    load(fname)
    index.res <- index.res + 1
    result.beta_m <- c(result.beta_m, beta_m)
  }
}

sd(result.beta_m)

l = quantile(result.beta_m, c(0.025, 0.975))

save(list = 'l', file = 'RData/summary.beta_m.RData')

