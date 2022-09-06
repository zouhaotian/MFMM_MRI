#setwd('/pine/scr/h/a/haotian/rev_code/')
options(mc.cores = parallel::detectCores())
library(tidyverse)
library(data.table)
library(splines)
library(survival)
library(Matrix)
library(rstan)
library(lme4)
library(mgcv)
library(sn)
library(tdROC)
library(ipred)
library(Amelia)
library(MBA)
library(RSpectra)
library(mfaces)
library(orthogonalsplinebasis)
library(statmod)
rstan_options(auto_write = FALSE)


model = 's1_MJM_MRI'

source('source_code/functions_MJM_MRI.R')
source('source_code/data_gen1.R')
load('RData/data_genx.RData')

# seed <- 1
N.train <- tp$N.train; N.test <- tp$N.test; N <- N.train + N.test
ID <- 1:N; ID.train <- 1:N.train; ID.test <- (N.train+1):N

dat <- sim(seed*2016, tp, N)
long <- dat$long; long.train <- long[long$ID %in% ID.train, ]; long.test <- long[long$ID %in% ID.test, ]
surv <- dat$surv; surv.train <- surv[surv$ID %in% ID.train, ]; surv.test <- surv[surv$ID %in% ID.test, ]
MRI.dat = dat$MRI.dat; M.train <- MRI.dat[ID.train, ]; M.test <- MRI.dat[ID.test, ]

J = ncol(long.train) - 2
g_obs = M.train
V = ncol(M.train)
obsgrid <- c(1:V)/V

## FPCA, and then create B matrices ##
FPCA.result <- FPCA(g_obs, obsgrid)
Lb <- ncol(FPCA.result$phi_est)
B1 <- create_B(g_obs, FPCA.result$mu_est, FPCA.result$phi_est, obsgrid, L=Lb, by=0.005)
B2 <- B1
B3 <- B1
B4 <- B1
B1.test <- create_B(M.test, FPCA.result$mu_est, FPCA.result$phi_est, obsgrid, L=Lb, by=0.005)
B2.test <- B1.test
B3.test <- B1.test
B4.test <- B1.test

nobs <- nrow(long.train)

##
tau <- c(0, 1)
n.tau <- length(tau) - 1

## get time spline matrix ##
Y.spline = matrix(NA, nrow = nobs, ncol = n.tau)
for (i in 1:nobs){
  for (k in 1:n.tau){
    Y.spline[i, k] = ifelse(long.train[i, 2]>tau[k], long.train[i, 2]-tau[k], 0)
  }
}

## get survival time spline matrix ##
surv.spline <- vector('list', 1)
tmp.time <- surv.train$surv_time
tmp.spline <- matrix(NA, length(tmp.time), n.tau)
for (i in 1:length(tmp.time)){
  tmp.spline[i, ] <- get_tau(tmp.time[i])
}
surv.spline[[1]] <- tmp.spline


## create h matrix (baseline hazard index matrix) ##
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

# get initial values #
# l1 <- calc_init_long1(long.train$ID, long.train[, 3], Y.spline, long.train[, 2], B1)
# l2 <- calc_init_long2(long.train$ID, long.train[, 4], Y.spline, B2, re = l1[[2]])
# l3 <- calc_init_long2(long.train$ID, long.train[, 5], Y.spline, B3, re = l1[[2]])
# rho <- l1$corr.est
# paras.est1 <- l1[[1]]
# paras.est2 <- l2[[1]]
# paras.est3 <- l3[[1]]
# beta1 <- paras.est1[1:(1+n.tau)]
# beta2 <- paras.est2[1:(1+n.tau)]
# beta3 <- paras.est3[1:(1+n.tau)]
# B1_est <- paras.est1[(2+n.tau):(1+n.tau+Lb)]
# B2_est <- paras.est2[(2+n.tau):(1+n.tau+Lb)]
# B3_est <- paras.est3[(2+n.tau):(1+n.tau+Lb)]
# v2_hat <- paras.est2[(2+n.tau+Lb):length(paras.est2)]
# v3_hat <- paras.est3[(2+n.tau+Lb):length(paras.est3)]
w <- as.matrix(surv.train$x)
beta1 <-  beta2 <- beta3 <- rep(0, 1+n.tau)
B1_est <- B2_est <- B3_est <- rep(0, Lb)
v2_hat <- v3_hat <- rep(0, 2)

## Determine missingness in Y ##
Y_missing <- missing_determine(value = list(Y1 = long.train[, 3],
                                            Y2 = long.train[, 4],
                                            Y3 = long.train[, 5]))
mean_Y <- list(mean1 = rep(mean(long.train[, 3], na.rm = T), Y_missing$len.missing[1]), 
               mean2 = rep(mean(long.train[, 4], na.rm = T), Y_missing$len.missing[2]),
               mean3 = rep(mean(long.train[, 5], na.rm = T), Y_missing$len.missing[3]))
Y1_imp = mean_Y[[1]]
Y2_imp = mean_Y[[2]]
Y3_imp = mean_Y[[3]]

G <- 15 ## number of Gaussian quadrature points
weights <- gauss.quad(G, kind = 'legendre')$weights ## w_g^* in document
nodes <- gauss.quad(G, kind = 'legendre')$nodes  ## x_g^* in document

##  initial values ##
md <- stan_model(file = 'source_code/MFJM_missing_rho.stan')
stanDat <- list(n = N.train, J = J, nobs = nobs, nmiss = Y_missing$len.missing,
                ID1 = long.train$ID, ID2 = long.train$ID, ID3 = long.train$ID,
                Lb = Lb, P_surv = ncol(w),
                Ltau = n.tau, tau = as.array(tau[1:n.tau]), 
                time_spline1 = Y.spline, time_spline2 = Y.spline, time_spline3 = Y.spline,
                B1 = B1, B2 = B2, B3 = B3, Bw = B4, 
                time1 = long.train$Time, 
                time2 = long.train$Time, 
                time3 = long.train$Time,
                w = w, 
                Y1 = Y_missing$new.value[[1]], 
                Y2 = Y_missing$new.value[[2]], 
                Y3 = Y_missing$new.value[[3]],
                miss_index1 = Y_missing$missing.index[[1]], 
                miss_index2 = Y_missing$missing.index[[2]], 
                miss_index3 = Y_missing$missing.index[[3]], 
                surv_time = surv.train$surv_time, status = surv.train$status,
                h_grid = h.grid, h_index = h.index, 
                zero = rep(0, 2))
pars <- c('beta1', 'beta2', 'beta3', 
          'sigma_u', 'rho', 'v2', 'v3', 'omega', 
          'logh0', 'gamma', 'alpha', 'BX1', 'BX2', 'BX3', 'BW',
          'Y1_imp', 'Y2_imp', 'Y3_imp', 'u')

rnd = 1/3
inits1 <- list(beta1=beta1+rnd, beta2=beta2+rnd, beta3=beta3+rnd,
               BX1=B1_est+rnd, BX2=B2_est+rnd, BX3=B3_est+rnd, BW = rep(0, Lb),
               sigma_u = rep(1, 2)+rnd, 
               rho = 0+rnd,
               v2=v2_hat+rnd, v3=v3_hat+rnd, 
               omega = rep(1, J)+rnd, 
               logh0 = as.array(rep(0, n.tau)), 
               gamma = as.array(0+rnd), alpha = rep(0, J)+rnd,
               Y1_imp = Y1_imp+rnd, Y2_imp = Y2_imp+rnd, Y3_imp = Y3_imp+rnd)
inits2 <- list(beta1=beta1-rnd, beta2=beta2-rnd, beta3=beta3-rnd,
               BX1=B1_est-rnd, BX2=B2_est-rnd, BX3=B3_est-rnd, BW = rep(0, Lb),
               sigma_u = rep(1, 2)-rnd, 
               rho = 0-rnd,
               v2=v2_hat-rnd, v3=v3_hat-rnd, 
               omega = rep(1, J)-rnd, 
               logh0 = as.array(rep(0, n.tau)), 
               gamma = as.array(0-rnd), alpha = rep(0, J)-rnd,
               Y1_imp = Y1_imp-rnd, Y2_imp = Y2_imp-rnd, Y3_imp = Y3_imp-rnd)

inits <- list(c1=inits1, c2=inits2)
fitStan <- sampling(md, data = stanDat, iter = 2000, warmup = 1000, 
                    chains = 2, thin=1, init = inits, pars = pars, seed = 123,
                    control = list(adapt_delta = 0.8, max_treedepth=10))
fname <- paste0('stan_fit/', model, '/', seed, '_posterior.RData')
save(list = 'fitStan', file=fname)

## Extract samples ##
s <- list()
s_mean <- list()
for (i in 1:15){
  s[[i]] <- extract(fitStan, par = pars[i])[[1]]
  names(s)[i] <- pars[i]
  if (length(dim(s[[i]]))==1) s_mean[[i]] <- mean(s[[i]]) 
  else s_mean[[i]] <- apply(s[[i]], 2, mean)
  names(s_mean)[i] <- pars[i]
}

## compute true survival probability at time t: exp(-cumulative hazard)
G <- 15 ## number of Gaussian quadrature points
weights <- gauss.quad(G, kind = 'legendre')$weights ## w_g^* in document
nodes <- gauss.quad(G, kind = 'legendre')$nodes  ## x_g^* in document
H_true <- function(tp, x, weights, nodes, 
                   xi, zeta1, zeta2, zeta3, xi_m, t){
  f <- function(t) exp(tp$gamma0*(xi[1]*sqrt(2)*sin(pi*t) + xi[2]*sqrt(2)*cos(pi*t)) + 
                         tp$gamma11*zeta1*sqrt(2)*cos(2*pi*t) + 
                         tp$gamma12*zeta2*sqrt(2)*cos(2*pi*t) + 
                         tp$gamma13*zeta3*sqrt(2)*cos(2*pi*t) )
  sum.approx <- sum(weights*f(t/2*nodes + t/2))
  H <- exp(tp$logh0 + x*tp$gamma_x + sum(tp$gamma_m*xi_m))*t/2*sum.approx
  return(H)
}

## Dynamic prediction ##
starting.time <- c(0.3, 0.4, 0.5, 0.55, 0.6) 
delta.time <- seq(0.1, 0.25, by = 0.01)
time.index <- 0
AUC <- BS <- matrix(NA, length(starting.time), length(delta.time))
AUC.true <- BS.true <- matrix(NA, length(starting.time), length(delta.time))

md_dp <- stan_model('source_code/DP_rho.stan')
pars_dp <- c('u', 'S', 'S2', 'cond_S')

for (Tstart in starting.time){
  ## Filter the survival testing data with subject survives up to Tstart
  ## Also filter the corresponding longitudinal data
  filter.obj <- filter.time(long.test, surv.test, Tstart)
  surv.filtered <- filter.obj$surv.filtered
  surv.filtered$newID <- 1:nrow(surv.filtered)
  
  long.filtered <- filter.obj$long.filtered
  long.filtered$newID <- rep(1:nrow(surv.filtered), table(long.filtered$ID))
  long.filtered2 <- long.filtered[which(long.filtered$Time<=Tstart), ]
  tmp.ID.test <- unique(long.filtered2$newID)
  nobs.test <- length(long.filtered2[, 3])
  
  ## Get true xi and xi_m ##
  true_xi <- dat$xi[surv.filtered$ID, ]
  true_zeta1 <- as.matrix(dat$zeta[[1]][surv.filtered$ID, ])
  true_zeta2 <- as.matrix(dat$zeta[[2]][surv.filtered$ID, ])
  true_zeta3 <- as.matrix(dat$zeta[[3]][surv.filtered$ID, ])
  true_xi_m <- dat$xi_m[surv.filtered$ID, ]
  
  Y_missing_test <- missing_determine(value = list(Y1 = long.filtered2[, 3],
                                                   Y2 = long.filtered2[, 4],
                                                   Y3 = long.filtered2[, 5]))
  tmp.B1 <- B1.test[surv.filtered$ID - N.train, ]
  tmp.B2 <- B2.test[surv.filtered$ID - N.train, ]
  tmp.B3 <- B3.test[surv.filtered$ID - N.train, ]
  tmp.B4 <- B4.test[surv.filtered$ID - N.train, ]
  
  ## get time spline matrix ##
  Y.spline.test = matrix(NA, nrow = nobs.test, ncol = n.tau)
  for (i in 1:nobs.test){
    for (k in 1:n.tau){
      Y.spline.test[i, k] = ifelse(long.filtered2[i, 2]>tau[k], long.filtered2[i, 2]-tau[k], 0)
    }
  }
  
  w.test <- as.matrix(surv.filtered$x)
  x.test = w.test
  
  for (Tdelta in delta.time){
    Tstop <- Tstart + Tdelta
    time.index <- time.index + 1
    
    ## 
    h.index.test <- h.index.test2 <- matrix(0, nrow = nrow(surv.filtered), ncol = n.tau)
    h.grid.test <- h.grid.test2 <- matrix(0, nrow = nrow(surv.filtered), ncol = n.tau)
    
    for (i in 1:nrow(surv.filtered)){
      tmp.survtime <- Tstart
      for (k in 1:n.tau){
        h.index.test[i, k] <- ifelse(tmp.survtime>=tau[k] & tmp.survtime<tau[k+1], 1, 0)
        h.grid.test[i, k] <- ifelse(tmp.survtime>=tau[k], min(tmp.survtime, tau[k+1]), tmp.survtime)
      }
      tmp.survtime <- Tstop
      for (k in 1:n.tau){
        h.index.test2[i, k] <- ifelse(tmp.survtime>=tau[k] & tmp.survtime<tau[k+1], 1, 0)
        h.grid.test2[i, k] <- ifelse(tmp.survtime>=tau[k], min(tmp.survtime, tau[k+1]), tmp.survtime)
      }
    }
    h.grid.test <- cbind(0, h.grid.test)
    h.grid.test2 <- cbind(0, h.grid.test2)
    
    ## Sample random effects ##
    stanDat <- list(n = nrow(surv.filtered), J = J, nobs = nobs.test, 
                    ID1 = long.filtered2$newID, 
                    ID2 = long.filtered2$newID, 
                    ID3 = long.filtered2$newID,
                    Lb = Lb, P_surv = ncol(w.test),
                    Ltau = n.tau, tau = as.array(tau[1:n.tau]), 
                    time_spline1 = Y.spline.test, 
                    time_spline2 = Y.spline.test, 
                    time_spline3 = Y.spline.test,
                    B1 = tmp.B1, B2 = tmp.B2, B3 = tmp.B3, Bw = tmp.B4, 
                    time1 = long.filtered2$Time, 
                    time2 = long.filtered2$Time, 
                    time3 = long.filtered2$Time,
                    w = w.test, 
                    Y1 = Y_missing_test$new.value[[1]], 
                    Y2 = Y_missing_test$new.value[[2]], 
                    Y3 = Y_missing_test$new.value[[3]],
                    surv_time = rep(Tstart, nrow(surv.filtered)), 
                    surv_time2 = rep(Tstop, nrow(surv.filtered)), 
                    status = rep(0, nrow(surv.filtered)),
                    h_grid = h.grid.test, h_index = h.index.test, 
                    h_grid2 = h.grid.test2, h_index2 = h.index.test2, 
                    zero = rep(0, 2), 
                    beta1 = s_mean$beta1, beta2 = s_mean$beta2, beta3 = s_mean$beta3, 
                    BX1 = s_mean$BX1, BX2 = s_mean$BX2, BX3 = s_mean$BX3, BW = s_mean$BW, 
                    sigma_u = s_mean$sigma_u, rho = s_mean$rho, 
                    v2 = s_mean$v2, v3 = s_mean$v3, 
                    omega = s_mean$omega, 
                    logh0 = as.array(s_mean$logh0), 
                    gamma = as.array(s_mean$gamma), 
                    alpha = s_mean$alpha)
    fitStan_dp <- sampling(md_dp, data = stanDat, iter = 2000, warmup = 1000, 
                           chains = 2, thin=1, pars = pars_dp, seed = 123,
                           control = list(adapt_delta = 0.8, max_treedepth=10))
    surv.est <- colMeans(extract(fitStan_dp, par = 'cond_S')$cond_S)
    
    ## True survival probability 
    surv.true <- rep(NA, length(tmp.ID.test))
    for (i in 1:length(tmp.ID.test)){
      xi_i_true <- true_xi[i, ]
      zeta1_i_true <- true_zeta1[i, ]
      zeta2_i_true <- true_zeta2[i, ]
      zeta3_i_true <- true_zeta3[i, ]
      xi_m_i_true <- true_xi_m[i, ]
      
      Hi_true_t0 <- H_true(tp, x.test[i, ], weights, nodes, 
                           xi_i_true, zeta1_i_true, zeta2_i_true, zeta3_i_true, 
                           xi_m_i_true, t = Tstart)
      
      Hi_true_t1 <- H_true(tp, x.test[i, ], weights, nodes, 
                           xi_i_true, zeta1_i_true, zeta2_i_true, zeta3_i_true, 
                           xi_m_i_true, t = Tstart+Tdelta)
      
      surv.true[i] <- exp(-Hi_true_t1 + Hi_true_t0)
    }
    
    ## AUC and BS
    ROC.est <- tdROC(X = 1 - surv.est, Y = surv.filtered$surv_time, 
                     delta = surv.filtered$status, tau = Tstop,
                     span = 0.1, alpha = 0.05,
                     n.grid = 1000, cut.off = 0.5)
    AUC[which(Tstart==starting.time), which(Tdelta==delta.time)] <- ROC.est$AUC$value
    ROC.true <- tdROC(X = 1 - surv.true, Y = surv.filtered$surv_time, 
                      delta = surv.filtered$status, tau = Tstop,
                      span = 0.1, alpha = 0.05,
                      n.grid = 1000, cut.off = 0.5)
    AUC.true[which(Tstart==starting.time), which(Tdelta==delta.time)] <- ROC.true$AUC$value
    
    ## Calculate BS ##
    surv.obj <- Surv(surv.filtered$surv_time, surv.filtered$status)
    BS.true[which(Tstart==starting.time), which(Tdelta==delta.time)] <- sbrier(surv.obj, surv.true, btime = Tstop)
    BS[which(Tstart==starting.time), which(Tdelta==delta.time)] <- sbrier(surv.obj, surv.est, btime = Tstop)
  }
  
}

fname <- paste0('result/', model, '/', seed, '_dp.RData')
save(list = c('AUC', 'BS', 'AUC.true', 'BS.true'), file = fname)
