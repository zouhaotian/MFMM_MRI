library(tidyverse)
library(gridExtra)
library(latex2exp)
library(xtable)
library(data.table)

model = 's1_MJM_MRI'

J <- 3
P <- 9 

datadir <- paste0('result/', model, '/')
result.AUC <- result.AUC.true <- result.BS <-  result.BS.true <- list()
index.res <- 0
for (i in 1:120){
  fname2 <- paste0(datadir, i, '_dp.RData')
  if (file.exists(fname2)){
    load(fname2)
    index.res <- index.res + 1
    
    result.AUC[[index.res]] <- AUC
    result.BS[[index.res]] <- BS
    result.AUC.true[[index.res]] <- AUC.true
    result.BS.true[[index.res]] <- BS.true
    if (index.res==100) break;
  }
}


## iAUC is calculated via integration over \delta_t: (0.1, 0.11, ..., 0.25)
## Integration is calculated by Simpson's Rule:
## For 13 points: (0.1, ..., 0.22), we adopt Composite Simpson's Rule
## For 4 points: (0.22, 0.23 ,0.24, 0.25), we adopt 3/8 Simpson's Rule
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


# For the remaining 4 points #
composite.simpson.2 <- function(vec.remain, grid.points){
  h <- grid.points[2] - grid.points[1]
  int.approx <- 3/8*h*(vec.remain[1] + 3*vec.remain[2] + 3*vec.remain[3] + vec.remain[4])
  return(int.approx) 
}

starting.time <- c(0.3, 0.4, 0.5, 0.55, 0.6) 
delta.time <- seq(0.1, 0.25, by = 0.01)
l.time <- length(starting.time)
l.delta.time <- length(delta.time)
tg.1 <- delta.time[1:(l.delta.time-3)]
tg.2 <- delta.time[(l.delta.time-3):(l.delta.time)]

AUC.matrix <- BS.matrix <- matrix(0, l.time, l.delta.time)
AUC.true.matrix <- BS.true.matrix <- matrix(0, l.time, l.delta.time)
integrated_AUC_BS <- matrix(0, l.time, 4)
index1 <- 0
for (i in 1:index.res){
  tmp.AUC <- result.AUC[[i]]
  tmp.AUC.true <- result.AUC.true[[i]]
  tmp.BS <- result.BS[[i]]
  tmp.BS.true <- result.BS.true[[i]]
  
  if (sum(is.na(tmp.AUC))==0 & sum(is.na(tmp.AUC.true))==0){
    index1 <- index1 + 1
    AUC.matrix <- AUC.matrix + tmp.AUC
    BS.matrix <- BS.matrix + tmp.BS
    AUC.true.matrix <- AUC.true.matrix + tmp.AUC.true
    BS.true.matrix <- BS.true.matrix + tmp.BS.true
    
    for (j in 1:l.time){
      integrated_AUC_BS[j, 1] <- integrated_AUC_BS[j, 1] + 
        composite.simpson(tmp.AUC.true[j, 1:(l.delta.time-3)], tg.1) + 
        composite.simpson.2(tmp.AUC.true[j, (l.delta.time-3):(l.delta.time)], tg.2)
      integrated_AUC_BS[j, 2] <- integrated_AUC_BS[j, 2] + 
        composite.simpson(tmp.AUC[j, 1:(l.delta.time-3)], tg.1) + 
        composite.simpson.2(tmp.AUC[j, (l.delta.time-3):(l.delta.time)], tg.2)
      integrated_AUC_BS[j, 3] <- integrated_AUC_BS[j, 3] + 
        composite.simpson(tmp.BS.true[j, 1:(l.delta.time-3)], tg.1) + 
        composite.simpson.2(tmp.BS.true[j, (l.delta.time-3):(l.delta.time)], tg.2)
      integrated_AUC_BS[j, 4] <- integrated_AUC_BS[j, 4] + 
        composite.simpson(tmp.BS[j, 1:(l.delta.time-3)], tg.1) + 
        composite.simpson.2(tmp.BS[j, (l.delta.time-3):(l.delta.time)], tg.2)
    }
    
  }
  
}

AUC.matrix <- AUC.matrix/index1
AUC.true.matrix <- AUC.true.matrix/index1
BS.matrix <- BS.matrix/index1
BS.true.matrix <- BS.true.matrix/index1
integrated_AUC_BS <- integrated_AUC_BS/index1/(max(delta.time) - min(delta.time))

rownames(AUC.matrix) <- rownames(AUC.true.matrix) <- rownames(BS.matrix) <- rownames(BS.true.matrix) <- paste0('Landmark Time = ', starting.time)
colnames(AUC.matrix) <- colnames(AUC.true.matrix) <- colnames(BS.matrix) <- colnames(BS.true.matrix) <- paste0('Window = ', delta.time)
rownames(integrated_AUC_BS) <- paste0('T = ', starting.time)
colnames(integrated_AUC_BS) <- c('True iAUC', 'Estimated iAUC', 'True iBS', 'Estimated iBS')

xtable(integrated_AUC_BS, digits = c(0, rep(3, 4)))

fname = paste0("RData/integrated_AUC_BS_", model, ".csv")
write.table(integrated_AUC_BS, file = fname, append = F, row.names = TRUE, col.names = TRUE, sep = ',')

