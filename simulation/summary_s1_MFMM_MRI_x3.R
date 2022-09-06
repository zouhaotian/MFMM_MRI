library(tidyverse)
library(gridExtra)
library(latex2exp)
library(xtable)
library(data.table)

model = 's1_MFMM_MRI_x3'

J <- 3
P <- 9 

datadir <- paste0('result/', model, '/')
result.parameters <- result.mean <- result.phi <- result.psi <- 
  result.phi_m <- result.psi_m <- result.xi <- result.xi_m <- list()
result.AUC <- result.AUC.true <- result.BS <-  result.BS.true <- list()
index.res <- 0
for (i in 1:120){
  fname <- paste0(datadir, i, '_posterior.RData')
  fname2 <- paste0(datadir, i, '_dp.RData')
  if (file.exists(fname) & file.exists(fname2)){
    load(fname)
    load(fname2)
    index.res <- index.res + 1
    result.parameters[[index.res]] <- summary.parameters
    result.mean[[index.res]] <- summary.mean
    result.phi[[index.res]] <- phi_est
    result.psi[[index.res]] <- psi_est
    result.phi_m[[index.res]] <- Phi_m_est
    result.psi_m[[index.res]] <- Psi_m_est
    result.xi[[index.res]] <- summary.xi
    result.xi_m[[index.res]] <- summary.xi_m
    
    result.AUC[[index.res]] <- AUC
    result.BS[[index.res]] <- BS
    result.AUC.true[[index.res]] <- AUC.true
    result.BS.true[[index.res]] <- BS.true
    if (index.res==100) break;
  }
}

get_summary <- function(result, row.index, true.value){
  tmp.mean <- rep(NA, length(result))
  tmp.se <- rep(NA, length(result))
  tmp.cp <- 0
  for (i in 1:length(result)){
    tmp.mean[i] <- result[[i]][row.index, 1]
    tmp.se[i] <- result[[i]][row.index, 2]
    if (result[[i]][row.index, 3]<true.value & result[[i]][row.index, 4]>true.value){
      tmp.cp <- tmp.cp + 1
    }
  }
  v <- c(mean = mean(tmp.mean), 
         bias = mean(tmp.mean) - true.value, 
         sd = sd(tmp.mean), 
         se = sqrt(sum(tmp.se^2)/length(result)),
         cp = tmp.cp/length(result))
  return(v)
}

load('RData/data_genx3.RData')

L0 <- tp$L0
L1 <- tp$L1
Lm <- tp$Lm

true_value <- c(tp$d0, tp$d1, tp$d_m, tp$beta[2:J], tp$beta_m, 
                tp$omega, tp$omega_m, 
                tp$logh0, tp$gamma_x,  
                tp$gamma0, tp$gamma11, tp$gamma12, tp$gamma13, 
                tp$gamma_m)
summ.dat <- data.frame(true.values = 0, mean = 0, bias = 0, sd = 0, se = 0, cp = 0)
for (i in 1:length(true_value)){
  summ.dat[i, 2:6] <- get_summary(result.parameters, i, true_value[i])
}
summ.dat$true.values <- true_value

ii = length(true_value)
summ.dat[(ii+1):(ii+L0), ] <- 0
for (i in 1:index.res){
  tmp.xi = result.xi[[i]]
  for (l in 1:L0){
    summ.dat[ii+l, 1] <- summ.dat[ii+l, 1] +
      sum((tmp.xi[, l] - tmp.xi[, L0+l])^2)/nrow(tmp.xi)
  }
}

summ.dat[(ii+1):(ii+L0), 1] <- summ.dat[(ii+1):(ii+L0), 1]/index.res

ii = length(true_value) + L0
summ.dat[(ii+1):(ii+Lm), ] <- 0
for (i in 1:index.res){
  tmp.xi_m = result.xi_m[[i]]
  for (l in 1:Lm){
    summ.dat[ii+l, 1] <- summ.dat[ii+l, 1] +
      sum((tmp.xi_m[, l] - tmp.xi_m[, Lm+l])^2)/nrow(tmp.xi_m)
  }
}

summ.dat[(ii+1):(ii+Lm), 1] <- summ.dat[(ii+1):(ii+Lm), 1]/index.res

rownames(summ.dat) <- c('d01', 'd02', 'd11', 
                        'd_m1', 'd_m2', 'd_m3', 'd_m4', 'd_m5', 
                        'beta_2', 'beta_3', 
                        'beta_m',
                        'sigma_1', 'sigma_2', 'sigma_3', 'sigma_m',
                        'logh0', 'gamma_x', 'gamma0',
                        'gamma11', 'gamma12', 'gamma13', 
                        'gamma_m1', 'gamma_m2', 'gamma_m3',
                        'gamma_m4', 'gamma_m5', 
                        'xi_{i1}', 'xi_{i2}', 
                        'xi_{mi1}', 'xi_{mi2}', 'xi_{mi3}', 
                        'xi_{mi4}', 'xi_{mi5}')

## Summary of functions ##
tnew <- tp$tnew
tg <- length(tnew)
V <- tp$V
obsgrid <- (1:V)/V

mu <- tp$mu
mu_m <- function(v) 1 + sin(2*pi*v) + cos(2*pi*v)

mu_tnew <- mu
mu_m_obs = mu_m(obsgrid)


## get estimated function and 95% CI, and plot #
func_est <- function(result, time.grid, function.index, index, B, lab.x, lab.y, range.y){
  mse <- rep(0, length(result))
  B_est <- matrix(NA, length(result), length(time.grid))
  B_summary <- data.frame(mean = rep(0, length(time.grid)), 
                          lower = rep(0, length(time.grid)), 
                          upper = rep(0, length(time.grid)))
  for (i in 1:length(result)){
    B_est[i, ] <- result[[i]][index, function.index]
    mse[i] <- sum((B - B_est[i, ])^2)
  }
  for (s in 1:length(time.grid)){
    B_summary$mean[s] <- mean(B_est[, s])
    B_summary$lower[s] <- quantile(B_est[, s], 0.025)
    B_summary$upper[s] <- quantile(B_est[, s], 0.975)
  }
  B_summary$true <- B
  p <- ggplot(data = B_summary, aes(x = time.grid, y = true)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'grey70') + 
    geom_line(aes(y = true), linetype = 'solid', size = 0.7) + 
    geom_line(aes(y = mean), linetype = 'dashed', col = 'blue', size = 0.7) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill=NA, size=1),
          plot.margin = margin(0.5, 0.25, 0, 0.25, "cm")) +
    xlab(TeX(lab.x)) +
    ylab(TeX(lab.y)) + 
    ylim(range.y[1], range.y[2])
  l <- list(p = p, mse = mean(mse)/length(result))
  return(l)
}

## get estimated function and 95% CI, and plot #
func_est <- function(result, time.grid, function.index, index, B, lab.x, lab.y, range.y){
  mse <- rep(0, length(result))
  B_est <- matrix(NA, length(result), length(time.grid))
  B_summary <- data.frame(mean = rep(0, length(time.grid)), 
                          lower = rep(0, length(time.grid)), 
                          upper = rep(0, length(time.grid)))
  for (i in 1:length(result)){
    B_est[i, ] <- result[[i]][index, function.index]
  }
  for (s in 1:length(time.grid)){
    B_summary$mean[s] <- mean(B_est[, s])
    B_summary$lower[s] <- quantile(B_est[, s], 0.025)
    B_summary$upper[s] <- quantile(B_est[, s], 0.975)
  }
  mse <- mean((B - B_summary$mean)^2)
  B_summary$true <- B
  p <- ggplot(data = B_summary, aes(x = time.grid, y = true)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'grey70') + 
    geom_line(aes(y = true), linetype = 'solid', size = 0.7) + 
    geom_line(aes(y = mean), linetype = 'dashed', col = 'blue', size = 0.7) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill=NA, size=1),
          plot.margin = margin(0.5, 0.25, 0, 0.25, "cm")) +
    xlab(TeX(lab.x)) +
    ylab(TeX(lab.y)) + 
    ylim(range.y[1], range.y[2])
  l <- list(p = p, mse = mse)
  return(l)
}

l1 <- func_est(result.mean, tnew, 1, 1:tg, mu_tnew[[1]], lab.x = 'Time', lab.y = '$mu_1(t)$', c(7, 11))
l2 <- func_est(result.mean, tnew, 1, 1:tg+tg, mu_tnew[[2]], lab.x = 'Time', lab.y = '$mu_2(t)$', c(6, 8))
l3 <- func_est(result.mean, tnew, 1, 1:tg+tg*2, mu_tnew[[3]], lab.x = 'Time', lab.y = '$mu_3(t)$', c(1, 4))
lm <- func_est(result.mean, obsgrid, 1, 1:V+tg*3, mu_m_obs, lab.x = 'Voxel Location', lab.y ='$mu_m(t)$', c(-1, 3))

amse.mu1 <- l1[[2]]
amse.mu2 <- l2[[2]]
amse.mu3 <- l3[[2]]
amse.mum <- lm[[2]]
# 
# setEPS()
# postscript("plot/sim_mu1.eps", height = 8, width = 12)
# grid.arrange(l1[[1]], l2[[1]], l3[[1]], 
#              lm[[1]], nrow = 2)
# dev.off()

phi <- list(
  f1 = function(t) sqrt(2)*sin(pi*t),
  f2 = function(t) sqrt(2)*cos(pi*t)
)

phi_obs <- list(f1 = phi[[1]](tnew), f2 = phi[[2]](tnew))

psi <- list(
  f1 = function(t) sqrt(2)*cos(2*pi*t)
)

psi_obs <- list(f1 = psi[[1]](tnew))

phi1 <- func_est(result.phi, tnew, 1, 1:tg, phi_obs[[1]], lab.x = 'Time', lab.y = '$\\phi_1(t)$', c(-1, 2))
phi2 <- func_est(result.phi, tnew, 2, 1:tg, phi_obs[[2]], lab.x = 'Time', lab.y = '$\\phi_2(t)$', c(-2, 2))
amse.phi1 <- phi1$mse
amse.phi2 <- phi2$mse

psi1 <- func_est(result.psi, tnew, 1, 1:tg, psi_obs[[1]], lab.x = 'Time', lab.y = '$\\psi_1(t)$', c(-2, 2))
amse.psi1 <- psi1$mse

## Eigenfunctions in m_i(v)
phi_m <- list(
  f1 = function(v) sqrt(2)*sin(pi*v/2),
  f2 = function(v) sqrt(2)*sin(3*pi*v/2)
)

psi_m <- list(
  f1 = function(v) sqrt(2)*cos(pi*v/2),
  f2 = function(v) sqrt(2)*cos(3*pi*v/2),
  f3 = function(v) sqrt(2)*cos(5*pi*v/2),
  f4 = function(v) sqrt(2)*cos(7*pi*v/2),
  f5 = function(v) sqrt(2)*cos(9*pi*v/2)
)

phi_m_obs <- matrix(NA, nrow = length(obsgrid), ncol = L0)
psi_m_obs <- matrix(NA, nrow = length(obsgrid), ncol = Lm)
for (l in 1:L0) phi_m_obs[, l] <- phi_m[[l]](obsgrid)
for (l in 1:Lm) psi_m_obs[, l] <- psi_m[[l]](obsgrid)

phi_m.obj <- amse.phi_m <- vector('list', L0)
psi_m.obj <- amse.psi_m <- vector('list', Lm)
for (l in 1:L0){
  lab.y = paste0('$\\phi_{m', l, '}(t)$') 
  phi_m.obj[[l]] <- func_est(result.phi_m, obsgrid, l, 1:V, phi_m_obs[, l], lab.x = 'Voxel', lab.y = lab.y, c(-3, 3))
  amse.phi_m[[l]] = phi_m.obj[[l]]$mse
}

for (l in 1:Lm){
  lab.y = paste0('$\\psi_{m', l, '}(t)$') 
  psi_m.obj[[l]] <- func_est(result.psi_m, obsgrid, l, 1:V, psi_m_obs[, l], lab.x = 'Voxel', lab.y = lab.y, c(-3, 3))
  amse.psi_m[[l]] = psi_m.obj[[l]]$mse
}

# 
# setEPS()
# postscript("plot/sim_phi_psi1.eps", height = 4)
# grid.arrange(phi1[[1]], phi2[[1]], psi1[[1]], nrow = 1)
# dev.off()
# 
# setEPS()
# postscript("plot/sim_phi_m_psi_m1.eps", height = 8, width = 10)
# grid.arrange(phi_m.obj[[1]]$p, phi_m.obj[[2]]$p, 
#              psi_m.obj[[1]]$p, psi_m.obj[[2]]$p, 
#              psi_m.obj[[3]]$p, psi_m.obj[[4]]$p, 
#              psi_m.obj[[5]]$p, 
#              nrow = 4)
# dev.off()


summ.dat[1:(J+1+L0+L1+L0+Lm), 7] <- c('AMSE.mu1', 'AMSE.mu2', 'AMSE.mu3', 'AMSE.mum', 
                                      'AMSE.phi1', 'AMSE.phi2', 'AMSE.psi1', 
                                      'AMSE.phi_m1', 'AMSE.phi_m2', 
                                      'AMSE.psi_m1', 'AMSE.psi_m2', 'AMSE.psi_m3', 'AMSE.psi_m4', 'AMSE.psi_m5')
summ.dat[1:(J+1+L0+L1+L0+Lm), 8] <- c(amse.mu1, amse.mu2, amse.mu3, amse.mum, 
                                      amse.phi1, amse.phi2, amse.psi1,
                                      amse.phi_m[[1]], amse.phi_m[[2]], 
                                      amse.psi_m[[1]], amse.psi_m[[2]],
                                      amse.psi_m[[3]], amse.psi_m[[4]],
                                      amse.psi_m[[5]])
fname = paste0("RData/summary_", model, ".csv")
write.csv(summ.dat, file = fname)
kk = summ.dat[c(1, 3, 4, 5, 6, 7, 8)]
kk[, 1] = paste0(rownames(kk), '=', round(kk[, 1], 2))
kk[, 6] = substr(kk[, 6], 6, str_length(kk[, 6]))
print(xtable(kk, type = "latex", digits = c(0, rep(3, 7))), include.rownames=FALSE)
kk[(L0+L1+1):(L0+L1+Lm), c(4, 5)] = ''
kk[(L0+L1+Lm+J), c(4, 5)] = ''


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

