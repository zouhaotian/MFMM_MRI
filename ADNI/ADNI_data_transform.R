library(tidyverse)
library(car)

long <- read.csv('dataset/ADNI_long.csv')
J <- 5
lambda <- rep(NA, J)
for (i in 1:J){
  tmp.y <- long[, i+2]
  tmp.y <- tmp.y[!is.na(tmp.y)]
  
  tmp.y.t <- boxCox(tmp.y ~ 1, family="yjPower", plotit = F)
  tmp.lambda <- tmp.y.t$x  
  tmp.ll <- tmp.y.t$y
  lambda[i] <- tmp.lambda[which.max(tmp.ll)]
}

lambda[4] <- 0

long.2 <- long

for (i in 1:J){
  tmp.y <- long.2[, i+2]
  if (lambda[i]==0){
    long.2[, i+2] <- log(tmp.y+1)
  } else {
    long.2[, i+2] <- ((tmp.y+1)^lambda[i] - 1)/lambda[i]
  }
}

long.2$y4 <- log(long$y4)

write.csv(long.2, file = 'dataset/ADNI_long_2.csv', row.names = F)
