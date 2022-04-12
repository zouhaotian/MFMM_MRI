library(tidyverse)
library(gridExtra)
library(latex2exp)
library(xtable)
library(data.table)

J <- 5
P <- 9 

datadir <- c('result2/M1_whole_dp/', 'result2/M1_dp/', 'result2/M1_MFMM_dp/', 
             'result2/M2_whole_dp/', 'result2/M2_dp/', 'result2/M2_MFMM_dp/')

result.AUC <-  result.BS <- vector('list', length = length(datadir))
index.res <- rep(0, length(datadir))

for (m in 1:length(datadir)){
  for (i in 1:120){
    fname <- paste0(datadir[m], i, '.RData')
    if (file.exists(fname)){
      load(fname)
      index.res[m] <- index.res[m] + 1
      result.AUC[[m]][[index.res[m]]] <- AUC
      result.BS[[m]][[index.res[m]]] <- BS
      if (index.res[m]==100) break;
    }
  }
}


## iAUC is calculated via integration over \delta_t
## Integration is calculated by Simpson's Rule:
## For 81 points: (0.2, 0.21, ..., 1), we adopt Composite Simpson's Rule
composite.simpson <- function(vec, grid.points){
  n <- length(grid.points) - 1
  h <- grid.points[2] - grid.points[1]
  int.approx <- 0
  
  ## vector for R starts with index 1. So for even index, multiply by 4.
  for (j in 1:(n/2)){
    int.approx <- int.approx + 4*vec[2*j]
  }
  for (j in 1:(n/2-1)){
    int.approx <- int.approx + 2*vec[2*j+1]
  }
  int.approx <- int.approx + vec[1] + vec[n+1]
  int.approx <- h/3*int.approx
  return(int.approx)
}


starting.time <- c(2, 2.5, 3, 3.5, 4) 
delta.time <- seq(0.5, 1.5, by = 0.01)
l.time <- length(starting.time)
l.delta.time <- length(delta.time)
tg <- seq(0.5, 1.5, by = 0.01)
time.index = 1:l.delta.time

AUC.matrix <- BS.matrix <- list(M1 = matrix(0, l.time, l.delta.time), 
                                M2 = matrix(0, l.time, l.delta.time),
                                M3 = matrix(0, l.time, l.delta.time),
                                M4 = matrix(0, l.time, l.delta.time), 
                                M5 = matrix(0, l.time, l.delta.time), 
                                M6 = matrix(0, l.time, l.delta.time))
integrated_AUC_BS <- matrix(0, l.time, 2*length(datadir))
for (m in 1:length(datadir)){
  index2 = 0
  for (i in 1:index.res[m]){
    tmp.AUC <- result.AUC[[m]][[i]]
    tmp.BS <- result.BS[[m]][[i]]
    
    if (sum(is.na(tmp.AUC))==0){
      index2 <- index2 + 1
      AUC.matrix[[m]] <- AUC.matrix[[m]] + tmp.AUC
      BS.matrix[[m]] <- BS.matrix[[m]] + tmp.BS
      
      for (j in 1:l.time){
        integrated_AUC_BS[j, m*2-1] <- integrated_AUC_BS[j, m*2-1] + 
          composite.simpson(tmp.AUC[j, time.index], tg)
        integrated_AUC_BS[j, m*2] <- integrated_AUC_BS[j, m*2] + 
          composite.simpson(tmp.BS[j, time.index], tg)
      }
    }
    
    
    
  }
  AUC.matrix[[m]] <- AUC.matrix[[m]]/index2
  BS.matrix[[m]] <- BS.matrix[[m]]/index2
  integrated_AUC_BS[, (m*2-1):(m*2)] <- integrated_AUC_BS[, (m*2-1):(m*2)]/index2/(max(tg) - min(tg))
  
  rownames(AUC.matrix[[m]]) <- rownames(BS.matrix[[m]]) <- paste0('Landmark Time = ', starting.time)
  colnames(AUC.matrix[[m]]) <- colnames(BS.matrix[[m]]) <- paste0('Window = ', delta.time)
  
  fname = paste0('result2/AUC_BS_M', m, '.csv')
  write.table(AUC.matrix[[m]], file = fname, append = F, row.names = TRUE, col.names = FALSE, sep = ',')
  write.table(BS.matrix[[m]], file = fname, append = T, row.names = TRUE, col.names = FALSE, sep = ',')
  
}


rownames(integrated_AUC_BS) <- paste0('T = ', starting.time)
colnames(integrated_AUC_BS) <- rep(c('Estimated iAUC', 'Estimated iBS'), length(datadir))

xtable(integrated_AUC_BS, digits = c(0, rep(3, ncol(integrated_AUC_BS))))


write.table(integrated_AUC_BS, file = 'result2/iAUC_iBS.csv', append = F, row.names = TRUE, col.names = TRUE, sep = ',')
