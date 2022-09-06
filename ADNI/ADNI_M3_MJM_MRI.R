## Computation of DIC, EAIC, EBIC, LOOIC, WAIC for M3_MJM_MRI ##
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
library(mgcv)
rstan_options(auto_write = FALSE)

set.seed(2020)

source('source_code/functions.R')
source('source_code/f2.R')

long <- read.csv('dataset/ADNI_long_2.csv') 
surv <- read.csv('dataset/ADNI_surv.csv')
MRI.dat = read.csv('dataset/Surv.mri.dat.csv') %>% as.matrix()

N <- nrow(surv)
J <- 5 ## number of longitudinal outcomes

## FPCA2 wraps FPCA function
FPCA2 <- function(dat){
  S <- ncol(dat)
  obsgrid <- 0:(S-1)
  by <- 1
  FPCA.obj <- FPCA(dat, obsgrid, pve = 0.99, by = by, L = 8)
  return(FPCA.obj)
}
Y1.FPCA <- FPCA2(MRI.dat)
save(list = c('Y1.FPCA'), file = 'RData/FPCA_result_MJM_MRI.RData')

## The minimum FPCA eigenvalues is 6, so we set Lb = 6 ##
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

B1 <- create_B2(MRI.dat, Y1.FPCA)
B2 <- B1
B3 <- B1
B4 <- B1

## Refine the time grid ##
tnew <- (0:101)/10
tg <- length(tnew)

## Impute the missing longitudinal outcomes, and average imputed datasets ##
M <- 5
long.impute <- amelia(x = long, m = M, idvars = 'ID', ts = 'Time', splinetime = 6)
long.impute.dataset <- long.impute$imputations
long.i <- long
for (j in 1:J){
  tmp.y <- rep(0, nrow(long))
  for (m in 1:M){
    tmp.y <- tmp.y + long.impute.dataset[[m]][, j+2]
  }
  long.i[, j+2] <- tmp.y/M
}

Y.list.o <- list(Y1 = long$y1, Y2 = long$y2, Y3 = long$y3, Y4 = long$y4, Y5 = long$y5)
Y.list.i <- list(Y1 = long.i$y1, Y2 = long.i$y2, Y3 = long.i$y3, Y4 = long.i$y4, Y5 = long.i$y5)

## get initial values of linear mixed model and survival model ##
x_surv <- as.matrix(surv[, 4:7])
l1 <- calc_init_long1(long$ID, long.i$y1, long$Time, long$Time, B = B1)
l2 <- calc_init_long2(long$ID, long.i$y2, long$Time, re = l1[[2]], B = B1)
l3 <- calc_init_long2(long$ID, long.i$y3, long$Time, re = l1[[2]], B = B1)
l4 <- calc_init_long2(long$ID, long.i$y4, long$Time, re = l1[[2]], B = B1)
l5 <- calc_init_long2(long$ID, long.i$y5, long$Time, re = l1[[2]], B = B1)
v <- list(v2 = l2$fixed[3:4], v3 = l3$fixed[3:4], v4 = l4$fixed[3:4], v5 = l5$fixed[3:4])
sigma_hat <- sqrt(c(l1$var.est[3], l2$var.est, l3$var.est, l4$var.est, l5$var.est))

fit.cox <- coxph(Surv(surv_time, status) ~ Age + Gender + Education + APOE4,
                 data = surv, method = 'breslow')
gamma_x_hat <- unname(fit.cox$coefficients)
baseline.cumulative.hazard <- basehaz(fit.cox)

tau <- c(0, 4, 6, 8, max(surv$surv_time)+1)
n.tau <- length(tau) - 1
h.index <- matrix(0, nrow = N, ncol = n.tau)
h.grid <- matrix(0, nrow = N, ncol = n.tau)

for (i in 1:N){
  tmp.survtime <- surv$surv_time[i]
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
mean_Y <- list(mean1 = rep(mean(long$y1, na.rm = T), Y_missing$len.missing[1]), 
               mean2 = rep(mean(long$y2, na.rm = T), Y_missing$len.missing[2]),
               mean3 = rep(mean(long$y3, na.rm = T), Y_missing$len.missing[3]), 
               mean4 = rep(mean(long$y4, na.rm = T), Y_missing$len.missing[4]), 
               mean5 = rep(mean(long$y5, na.rm = T), Y_missing$len.missing[5]))
Y.list.new <- miss.index <- list()
for (j in 1:J){
  Y.list.new[[j]] <- Y_missing$new.value[[j]]
  miss.index[[j]] <- Y_missing$missing.index[[j]]
}

start_end_index <- index_determine(long$ID)

## Stan fit ##
md = stan_model('source_code/M3_MJM_MRI.stan')
stan_dat <- list(n = N, J = J, 
                 nobs = nrow(long), nmiss = Y_missing$len.missing, 
                 id_long = long$ID, 
                 P_surv = ncol(x_surv), Lb = Lb, 
                 time = long$Time, B1 = B1, 
                 w = x_surv, 
                 Y1 = Y.list.new[[1]], Y2 = Y.list.new[[2]], Y3 = Y.list.new[[3]],
                 Y4 = Y.list.new[[4]], Y5 = Y.list.new[[5]],
                 miss_index1 = miss.index[[1]], miss_index2 = miss.index[[2]], 
                 miss_index3 = miss.index[[3]], miss_index4 = miss.index[[4]], 
                 miss_index5 = miss.index[[5]],
                 surv_time = surv$surv_time, status = surv$status,
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
          'logh0', 'gamma', 'alpha',
          'u', 'LL', 'log_lik',
          'Y1_imp', 'Y2_imp', 'Y3_imp', 'Y4_imp', 'Y5_imp')
fitStan <- sampling(md, data = stan_dat, iter = 3000, warmup = 2000, 
                    chains = 2, thin=1, init = inits, pars = pars, seed = 2020,
                    control = list(adapt_delta = 0.8, max_treedepth=10))
fname <- paste0('stan_fit/', 'M3_MJM_MRI_posterior.RData')
save(list = 'fitStan', file = fname)
# load(fname)
summ <- summary(fitStan)$summary[1:(J*2 + Lb*6 + 2 + 1 + 2*(J-1) + J + n.tau + ncol(x_surv) + J), c(1, 3, 4, 8, 9, 10)]
write.csv(summ, file = 'RData/summary_posterior_M3_MJM_MRI.csv')

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
u_sample <- extract(fitStan, pars = 'u')$u ## Q*L0*N
LL_sample <- extract(fitStan, pars = 'LL')$LL
log_lik_sample <- extract(fitStan, pars = 'log_lik')$log_lik
Y_imp_sample <- extract(fitStan, pars = c('Y1_imp', 'Y2_imp', 'Y3_imp', 'Y4_imp', 'Y5_imp'))

ll_surv <- LL_sample

## Calculation of log likelihood ##
ll_long_full <- log_lik_sample

Dbar.long <- sum(-2*ll_long_full)/Q
Dbar.surv <- sum(-2*ll_surv)/Q
Dbar <- Dbar.long + Dbar.surv

np <- nrow(summ)
EAIC <- Dbar + 2*np
EBIC <- Dbar + log(N)*np

## Compute LOOIC, WAIC ##
ll_full <- ll_long_full + ll_surv
rel_n_eff <- relative_eff(exp(ll_full), chain_id = rep(1:2, each = 1000))
looic <- loo(ll_full, r_eff = rel_n_eff, cores = 4)$estimates[3, 1]
waic <- waic(ll_full)$estimates[3, 1]

l <- list(Dbar.long = Dbar.long, Dbar.surv = Dbar.surv,
          np = np, EAIC = EAIC, EBIC = EBIC,
          looic = looic, waic = waic)
save(list = 'l', file = 'RData/summary_posterior_M3_MJM_MRI_MP.RData')

