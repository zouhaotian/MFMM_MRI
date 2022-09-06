## Dynamic Prediction for M3_MJM_MRI ##
# setwd('/ADNI_2/')
options(mc.cores = parallel::detectCores())
library(tidyverse)
library(data.table)
library(splines)
library(lme4)
library(survival)
library(rstan)
library(Matrix)
library(refund)
library(orthogonalsplinebasis)
library(mfaces)
library(loo)
library(Amelia)
library(sn)
library(dlm)
library(tdROC)
library(ipred)
library(mgcv)
rstan_options(auto_write = FALSE)

# seed <- 1
set.seed(2028*seed)

source('source_code/functions.R')
source('source_code/f2.R')

long <- read.csv('dataset/ADNI_long_2.csv') 
surv <- read.csv('dataset/ADNI_surv.csv')
Y1.mri.dat <- fread('dataset/Surv.mri.dat.csv') %>% as.matrix()

N <- nrow(surv)
J <- 5 ## number of longitudinal outcomes

N.train <- round(N*0.7); N.test <- N - N.train; ID <- 1:N
ID.train <- sample(ID, size = N.train) %>% sort()
ID.test <- ID[-ID.train]

long.train <- long[long$ID %in% ID.train, ]
long.train$ID <- rep(1:N.train, table(long.train$ID))
long.test <- long[long$ID %in% ID.test, ]
long.test$ID <- rep((N.train+1):N, table(long.test$ID))

surv.train <- surv[surv$ID %in% ID.train, ]
surv.train$ID <- 1:N.train
surv.test <- surv[surv$ID %in% ID.test, ]
surv.test$ID <- (N.train+1):N

##
FPCA2 <- function(dat){
  S <- ncol(dat)
  obsgrid <- 0:(S-1)
  by <- 1
  FPCA.obj <- FPCA(dat, obsgrid, pve = 0.99, by = by, L = 8)
  return(FPCA.obj)
}
Y1.FPCA <- FPCA2(Y1.mri.dat[ID.train, ])
Lb <- ncol(Y1.FPCA$phi_est)

create_B2 <- function(dat, FPCA.obj){
  eigen.values <- FPCA.obj$values
  S <- ncol(dat)
  obsgrid <- 0:(S-1)/eigen.values[1]
  by <- 1/eigen.values[1]
  mu_est <- FPCA.obj$mu_est
  phi_est <- FPCA.obj$phi_est/sqrt(by)
  eigen.values <- eigen.values/eigen.values[1]
  B <- create_B(dat, mu_est, phi_est, obsgrid, L = Lb, by = by)
  return(B)
}

B1 <- create_B2(Y1.mri.dat[ID.train, ], Y1.FPCA)
B2 = B1
B3 = B1
B4 = B1

B1.test <- create_B2(Y1.mri.dat[ID.test, ], Y1.FPCA)

## Refine the time grid ##
tnew <- (0:101)/10
tg <- length(tnew)

## Impute the missing longitudinal outcomes, and average imputed datasets ##
M <- 5
long.impute <- amelia(x = long.train, m = M, idvars = 'ID', ts = 'Time',  splinetime = 6)
long.impute.dataset <- long.impute$imputations
long.i <- long.train
for (j in 1:J){
  tmp.y <- rep(0, nrow(long.train))
  for (m in 1:M){
    tmp.y <- tmp.y + long.impute.dataset[[m]][, j+2]
  }
  long.i[, j+2] <- tmp.y/M
}

Y.list.o <- list(Y1 = long.train$y1, Y2 = long.train$y2, Y3 = long.train$y3, Y4 = long.train$y4, Y5 = long.train$y5)
Y.list.i <- list(Y1 = long.i$y1, Y2 = long.i$y2, Y3 = long.i$y3, Y4 = long.i$y4, Y5 = long.i$y5)

## get initial values of linear mixed model and survival model ##
x_surv <- as.matrix(surv.train[, 4:7])
l1 <- calc_init_long1(long.train$ID, long.i$y1, long.train$Time, long.train$Time, B = B1)
l2 <- calc_init_long2(long.train$ID, long.i$y2, long.train$Time, re = l1[[2]], B = B1)
l3 <- calc_init_long2(long.train$ID, long.i$y3, long.train$Time, re = l1[[2]], B = B1)
l4 <- calc_init_long2(long.train$ID, long.i$y4, long.train$Time, re = l1[[2]], B = B1)
l5 <- calc_init_long2(long.train$ID, long.i$y5, long.train$Time, re = l1[[2]], B = B1)
v <- list(v2 = l2$fixed[3:4], v3 = l3$fixed[3:4], v4 = l4$fixed[3:4], v5 = l5$fixed[3:4])
sigma_hat <- sqrt(c(l1$var.est[3], l2$var.est, l3$var.est, l4$var.est, l5$var.est))

fit.cox <- coxph(Surv(surv_time, status) ~ Age + Gender + Education + APOE4,
                 data = surv.train, method = 'breslow')
gamma_x_hat <- unname(fit.cox$coefficients)
baseline.cumulative.hazard <- basehaz(fit.cox)

tau <- c(0, 4, 6, 8, max(surv.train$surv_time)+1)
n.tau <- length(tau) - 1
h.index <- matrix(0, nrow = N.train, ncol = n.tau)
h.grid <- matrix(0, nrow = N.train, ncol = n.tau)

for (i in 1:N.train){
  tmp.survtime <- surv.train$surv_time[i]
  for (k in 1:n.tau){
    h.index[i, k] <- ifelse(tmp.survtime>=tau[k] & tmp.survtime<tau[k+1], 1, 0)
    h.grid[i, k] <- ifelse(tmp.survtime>=tau[k], min(tmp.survtime, tau[k+1]), tmp.survtime)
  }
}
h.grid <- cbind(0, h.grid)

logh0_hat <- rep(0, n.tau)
for (k in 1:n.tau){
  index1 <- min(which(baseline.cumulative.hazard$time>=tau[k]))
  index2 <- max(which(baseline.cumulative.hazard$time<tau[k+1]))
  lm.obj2 <- lm(hazard ~ time, data = baseline.cumulative.hazard[index1:index2, ])
  logh0_hat[k] <- log(unname(lm.obj2$coefficients[2]))
}

## Determine missingness
Y_missing <- missing_determine(value = Y.list.o)
mean_Y <- list(mean1 = rep(mean(long.train$y1, na.rm = T), Y_missing$len.missing[1]), 
               mean2 = rep(mean(long.train$y2, na.rm = T), Y_missing$len.missing[2]),
               mean3 = rep(mean(long.train$y3, na.rm = T), Y_missing$len.missing[3]), 
               mean4 = rep(mean(long.train$y4, na.rm = T), Y_missing$len.missing[4]), 
               mean5 = rep(mean(long.train$y5, na.rm = T), Y_missing$len.missing[5]))
Y.list.new <- miss.index <- list()
for (j in 1:J){
  Y.list.new[[j]] <- Y_missing$new.value[[j]]
  miss.index[[j]] <- Y_missing$missing.index[[j]]
}

start_end_index <- index_determine(long.train$ID)

## Stan fit ##
md = stan_model('source_code/M3_MJM_MRI.stan')
stan_dat <- list(n = N.train, J = J, 
                 nobs = nrow(long.train), nmiss = Y_missing$len.missing, 
                 id_long = long.train$ID, Lb = Lb, 
                 P_surv = ncol(x_surv),
                 time = long.train$Time,  B1 = B1,
                 w = x_surv, 
                 Y1 = Y.list.new[[1]], Y2 = Y.list.new[[2]], Y3 = Y.list.new[[3]],
                 Y4 = Y.list.new[[4]], Y5 = Y.list.new[[5]],
                 miss_index1 = miss.index[[1]], miss_index2 = miss.index[[2]], 
                 miss_index3 = miss.index[[3]], miss_index4 = miss.index[[4]], 
                 miss_index5 = miss.index[[5]],
                 surv_time = surv.train$surv_time, status = surv.train$status,
                 Ltau = n.tau, tau = tau[1:n.tau], 
                 h_grid = h.grid, h_index = h.index, 
                 zero = rep(0, 2),
                 s1 = start_end_index$start, e1 = start_end_index$end)

rnd <- 1/3
inits1 <- list(beta1 = rep(0, 2) + rnd, beta2 = rep(0, 2) + rnd, 
               beta3 = rep(0, 2) + rnd, beta4 = rep(0, 2) + rnd, 
               beta5 = rep(0, 2) + rnd, 
               BX1 = rep(0, Lb) + rnd, BX2 = rep(0, Lb) + rnd, 
               BX3 = rep(0, Lb) + rnd, BX4 = rep(0, Lb) + rnd, 
               BX5 = rep(0, Lb) + rnd, BW = rep(0, Lb) + rnd, 
               sigma_u = sqrt(l1$var.est[1:2]), 
               rho = l1$corr.est,
               v2 = rep(0, 2) + rnd, v3 = rep(0, 2) + rnd, 
               v4 = rep(0, 2) + rnd, v5 = rep(0, 2) + rnd,
               omega = g.inits.sd(sigma_hat),
               logh0 = g.inits(logh0_hat), 
               gamma = as.array(g.inits(gamma_x_hat)), 
               alpha = c(1, -1, -1, -1, 1),
               Y1_imp = mean_Y[[1]]+rnd, Y2_imp = mean_Y[[2]]+rnd, Y3_imp = mean_Y[[3]]+rnd, 
               Y4_imp = mean_Y[[4]]+rnd, Y5_imp = mean_Y[[5]]+rnd)
inits2 <- list(beta1 = rep(0, 2) - rnd, beta2 = rep(0, 2) - rnd, 
               beta3 = rep(0, 2) - rnd, beta4 = rep(0, 2) - rnd, 
               beta5 = rep(0, 2) - rnd, 
               BX1 = rep(0, Lb) - rnd, BX2 = rep(0, Lb) - rnd, 
               BX3 = rep(0, Lb) - rnd, BX4 = rep(0, Lb) - rnd, 
               BX5 = rep(0, Lb) - rnd, BW = rep(0, Lb) - rnd, 
               sigma_u = sqrt(l1$var.est[1:2]), 
               rho = l1$corr.est,
               v2 = rep(0, 2) - rnd, v3 = rep(0, 2) - rnd, 
               v4 = rep(0, 2) - rnd, v5 = rep(0, 2) - rnd,
               omega = g.inits.sd(sigma_hat),
               logh0 = g.inits(logh0_hat), 
               gamma = as.array(g.inits(gamma_x_hat)), 
               alpha = c(1, -1, -1, -1, 1),
               Y1_imp = mean_Y[[1]]-rnd, Y2_imp = mean_Y[[2]]-rnd, Y3_imp = mean_Y[[3]]-rnd, 
               Y4_imp = mean_Y[[4]]-rnd, Y5_imp = mean_Y[[5]]-rnd)

inits <- list(c1 = inits1, c2 = inits2)
pars <- c('beta1', 'beta2', 'beta3', 'beta4', 'beta5',
          'BX1', 'BX2', 'BX3', 'BX4', 'BX5', 'BW', 
          'sigma_u', 'rho', 'v2', 'v3', 'v4', 'v5',
          'omega', 
          'logh0', 'gamma', 'alpha')
n.iters = 3000
n.warmups = 2000

fitStan <- sampling(md, data = stan_dat, iter = n.iters, warmup = n.warmups, 
                    chains = 2, thin=1, init = inits, pars = pars, seed = 2021,
                    control = list(adapt_delta = 0.8, max_treedepth=10))
fname = paste0('stan_fit/M3_MJM_MRI_dp/', seed, '_posterior.RData')
save(list = 'fitStan', file = fname)


## Extract samples ##
Q <- 2000
beta_sample <- extract(fitStan, pars = c('beta1', 'beta2', 'beta3', 'beta4', 'beta5'))
B_sample <- extract(fitStan, pars = c('BX1', 'BX2', 'BX3', 'BX4', 'BX5', 'BW'))
sigma_u_sample <- extract(fitStan, par = 'sigma_u')$sigma_u
rho_sample <- extract(fitStan, par = 'rho')$rho
v_sample <- extract(fitStan, pars = c('v2', 'v3', 'v4', 'v5'))
omega_sample <- extract(fitStan, par = 'omega')$omega
logh0_sample <- extract(fitStan, 'logh0')$logh0
gamma_sample <- extract(fitStan, pars = 'gamma')$gamma
alpha_sample <- extract(fitStan, pars = 'alpha')$alpha

# 
beta <- lapply(1:J, FUN = function(x) colMeans(beta_sample[[x]]))
B <- lapply(1:(J+1), FUN = function(x) colMeans(B_sample[[x]]))
sigma_u <- apply(sigma_u_sample, 2, 'mean')
rho <- mean(rho_sample)
v <- lapply(1:(J-1), FUN = function(x) colMeans(v_sample[[x]]))
omega <- apply(omega_sample, 2, 'mean')
logh0 <- apply(logh0_sample, 2, 'mean')
gamma <- apply(gamma_sample, 2, 'mean')
alpha <- apply(alpha_sample, 2, 'mean')

## Dynamic prediction settings
md_dp = stan_model('source_code/M3_MJM_MRI_dp.stan')
pars_dp = c('u')
starting.time <- c(2, 2.5, 3, 3.5, 4)
delta.time <- seq(0.5, 1.5, by = 0.01)
min.log <- log(.Machine$double.xmin) + 10
AUC <- BS <- matrix(NA, length(starting.time), length(delta.time))

## compute estimated survival probability at time t: exp(-cumulative hazard)
H_est <- function(B1, w, h_grid, u, t){
  H <- rep(NA, n.tau)
  
  
  g1 = beta[[1]][1] + sum(B1*B[[1]]) + u[1];
  g2 = beta[[2]][1] + sum(B1*B[[2]]) + v[[1]][1]*u[1];
  g3 = beta[[3]][1] + sum(B1*B[[3]]) + v[[2]][1]*u[1];
  g4 = beta[[4]][1] + sum(B1*B[[4]]) + v[[3]][1]*u[1];
  g5 = beta[[5]][1] + sum(B1*B[[5]]) + v[[4]][1]*u[1];
  
  rs_1 = alpha[1]*g1;
  rs_2 = alpha[2]*g2;
  rs_3 = alpha[3]*g3;
  rs_4 = alpha[4]*g4;
  rs_5 = alpha[5]*g5;
  
  coef_t = alpha[1]*(beta[[1]][2]+u[2]) + alpha[2]*(beta[[2]][2]+v[[1]][2]*u[2]) + 
    alpha[3]*(beta[[3]][2]+v[[2]][2]*u[2]) + alpha[4]*(beta[[4]][2]+v[[3]][2]*u[2]) + 
    alpha[5]*(beta[[5]][2]+v[[4]][2]*u[2]);
  
  rs_w = sum(w*gamma)
  
  for (k in 1:n.tau){
    H[k] = exp(logh0[k] + rs_w + rs_1 + rs_2 + rs_3 + rs_4 + rs_5)*
      (exp(coef_t*h_grid[k+1]) - exp(coef_t*h_grid[k]))/coef_t;
  }
  
  return(sum(H))
}

for (Tstart in starting.time){
  
  ## filter out subjects with observed survival time greater than starting time
  surv.test2 <- surv.test[surv.test$surv_time>Tstart, ]
  long.test2 <- long.test[long.test$ID %in% surv.test2$ID, ]
  surv.test2$newID <- 1:nrow(surv.test2)
  long.test2$newID <- rep(1:nrow(surv.test2), table(long.test2$ID))
  x.test.surv <- as.matrix(surv.test2[, 4:7])
  
  ## filter out longitudinal observations prior or equal to starting time
  long.test.prior <- long.test2[long.test2$Time<=Tstart, ]
  long.test.posterior <- long.test2[long.test2$Time>Tstart, ]
  tmp.ID.test <- unique(long.test2$newID)
  
  ## Determine missingness
  Y.list.test <- list(y1 = long.test.prior$y1, y2 = long.test.prior$y2, y3 = long.test.prior$y3, 
                      y4 = long.test.prior$y4, y5 = long.test.prior$y5)
  Y_missing.test <- missing_determine(value = Y.list.test)
  
  tmp.B1 <- B1.test[surv.test2$ID - N.train, ]
  
  h_grid_T0 <- rep(NA, n.tau)
  for (k in 1:n.tau){
    h_grid_T0[k] <- ifelse(Tstart>=tau[k], min(Tstart, tau[k+1]), Tstart)
  }
  h_grid_T0 <- c(0, h_grid_T0)
  
  # NUTS
  stan_dat_dp <- list(n = length(tmp.ID.test), J = J, 
                      nobs = nrow(long.test.prior),
                      id_long = long.test.prior$newID,
                      P_surv = ncol(x.test.surv), Lb = Lb, 
                      time = long.test.prior$Time, B1 = tmp.B1, 
                      w = x.test.surv,
                      Y1 = Y_missing.test$new.value[[1]], 
                      Y2 = Y_missing.test$new.value[[2]], 
                      Y3 = Y_missing.test$new.value[[3]],
                      Y4 = Y_missing.test$new.value[[4]],
                      Y5 = Y_missing.test$new.value[[5]],
                      Ltau = n.tau, tau = tau[1:n.tau], 
                      h_grid = h_grid_T0, zero = rep(0, 2),
                      beta1 = beta[[1]], beta2 = beta[[2]], beta3 = beta[[3]], 
                      beta4 = beta[[4]], beta5 = beta[[5]],
                      BX1 = B[[1]], BX2 = B[[2]], BX3 = B[[3]], BX4 = B[[4]], 
                      BX5 = B[[5]], BW = B[[6]], 
                      sigma_u = sigma_u, rho = rho,
                      v2 = v[[1]], v3 = v[[2]], v4 = v[[3]], v5 = v[[4]], 
                      omega = unlist(omega),
                      logh0 = logh0, gamma = gamma, alpha = alpha)
  inits1_dp <- inits2_dp <- list(u = matrix(0, nrow = nrow(surv.test2), ncol = 2))
  inits_dp <- list(c1 = inits1_dp, c2 = inits2_dp)
  fitStan_dp <- sampling(md_dp, data = stan_dat_dp, iter = 2000, warmup = 1000, 
                         chains = 2, thin=1, init = inits_dp, pars = pars_dp, seed = 123,
                         control = list(adapt_delta = 0.8, max_treedepth=10))
  u_sample <- extract(fitStan_dp, pars = 'u')$u ## extract u
  
  ## calculate conditional survival probabilities at T+\delta_T
  for (Tdelta in delta.time){
    
    h_grid_T1 <- rep(NA, n.tau)
    for (k in 1:n.tau){
      h_grid_T1[k] <- ifelse(Tstart+Tdelta>=tau[k], min(Tstart+Tdelta, tau[k+1]), Tstart+Tdelta)
    }
    h_grid_T1 <- c(0, h_grid_T1)
    
    tmp.surv.predict <- rep(NA, length(tmp.ID.test))
    
    for (i in 1:length(tmp.ID.test)){
      u_i <- u_sample[, i, ] 
      Hi_est_t0 <- Hi_est_t1 <- rep(NA, Q)
      for (q in 1:Q){
        Hi_est_t0[q] <- H_est(tmp.B1[i, ], x.test.surv[i, ], 
                              h_grid_T0, u_i[q, ], 
                              t = Tstart)
        Hi_est_t1[q] <- H_est(tmp.B1[i, ], x.test.surv[i, ], 
                              h_grid_T1, u_i[q, ], 
                              t = Tstart + Tdelta)
      }
      
      cond_surv_prob <- exp(-Hi_est_t1 + Hi_est_t0)
      tmp.surv.predict[i] <- mean(cond_surv_prob)
    }
    
    ROC.est <- tdROC(X = 1 - tmp.surv.predict, Y = surv.test2$surv_time,
                     delta = surv.test2$status, tau = Tstart + Tdelta,
                     span = 0.1, alpha = 0.05,
                     n.grid = 1000, cut.off = 0.5)
    AUC[which(Tstart==starting.time), which(Tdelta==delta.time)] <- ROC.est$AUC$value
    
    surv.obj <- Surv(surv.test2$surv_time, surv.test2$status)
    BS[which(Tstart==starting.time), which(Tdelta==delta.time)] <- sbrier(surv.obj, tmp.surv.predict, btime = Tstart + Tdelta)
  }
}

fname <- paste0('result/M3_MJM_MRI_dp/', seed, '.RData')
save(list = c('AUC', 'BS'), file = fname)
