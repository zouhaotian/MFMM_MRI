# setwd('/sim/')
options(mc.cores = parallel::detectCores())
library(tidyverse)
library(data.table)
library(splines)
library(survival)
library(Matrix)
library(rstan)
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

model = 's1_NM'

source('source_code/functions.R')
source('source_code/data_gen1.R')
source('source_code/sFPCA.R')
load('RData/data_genx.RData')

# seed <- 1
N.train <- tp$N.train; N.test <- tp$N.test; N <- N.train + N.test
ID <- 1:N; ID.train <- 1:N.train; ID.test <- (N.train+1):N

dat <- sim(seed*2016, tp, N)
long <- dat$long; long.train <- long[long$ID %in% ID.train, ]; long.test <- long[long$ID %in% ID.test, ]
surv <- dat$surv; surv.train <- surv[surv$ID %in% ID.train, ]; surv.test <- surv[surv$ID %in% ID.test, ]
MRI.dat = dat$MRI.dat; M.train <- MRI.dat[ID.train, ]; M.test <- MRI.dat[ID.test, ]

J <- tp$J ## number of longitudinal outcomes
V <- ncol(MRI.dat)
obsgrid = (1:V)/V
P <- 9 ## number of basis functions for B-spline
L0 <- tp$L0 ## npc for U_i(t)
L1 <- tp$L1 ## npc for W_ij(t)
Lm <- tp$Lm
tnew <- tp$tnew
tg <- length(tnew)

## fit Cox model for survival data ##
fit.cox <- coxph(Surv(surv_time, status) ~ x, data = surv.train, method = 'breslow')
gamma_x_hat <- unname(fit.cox$coefficients)
baseline.cumulative.hazard <- basehaz(fit.cox)
lm.obj <- lm(hazard ~ time, data = baseline.cumulative.hazard)
logh0_hat <- log(unname(lm.obj$coefficients[2]))
x <- as.matrix(surv.train$x)
P_surv = ncol(x)

## Impute the missing longitudinal outcomes, and average imputed datasets ##
M <- 5
long.impute <- amelia(x = long.train, m = M, idvars = 'ID', ts = 'Time', splinetime = 6)
long.impute.dataset <- long.impute$imputations
long.i <- long.train
for (j in 1:J){
  tmp.y <- rep(0, nrow(long.train))
  for (m in 1:M){
    tmp.y <- tmp.y + long.impute.dataset[[m]][, j+2]
  }
  long.i[, j+2] <- tmp.y/M
}

Y.list.o <- list(Y1 = long.train$y1, Y2 = long.train$y2, Y3 = long.train$y3)
Y.list.i <- list(Y1 = long.i$y1, Y2 = long.i$y2, Y3 = long.i$y3)

dat.mface <- list('y1' = data.frame('subj' = long.train$ID, 'argvals' = long.train$Time, 'y' = Y.list.i[[1]]),
                  'y2' = data.frame('subj' = long.train$ID, 'argvals' = long.train$Time, 'y' = Y.list.i[[2]]),
                  'y3' = data.frame('subj' = long.train$ID, 'argvals' = long.train$Time, 'y' = Y.list.i[[3]]))
fit.mface <- mface.sparse(dat.mface, argvals.new = tnew, knots = 6, newdata = dat.mface)
# save(list = 'fit.mface', file = 'RData/fit_mface.RData')

## 
mu_est <- list()
long.demean <- long.i
for (j in 1:J){
  tmp.mu <- fit.mface$fit[[j]]$mu.new
  l <- bs.smooth(tmp.mu, tnew, tnew, nbasis = P)
  mu_est[[j]] <- list(value = tmp.mu, 
                      argvals = tnew, 
                      coefficient = l$coef.est)
  
  l <- bs.smooth(tmp.mu, tnew, long.train[, 2], nbasis = P)
  long.demean[, j+2] = long.i[, j+2] - l$est.value
}


C <- as.matrix(fit.mface$Chat.new)
C11 <- C[1:tg, 1:tg]
C12 <- C[1:tg, 1:tg+tg]
C13 <- C[1:tg, 1:tg+tg*2]
C22 <- C[1:tg+tg, 1:tg+tg]
C23 <- C[1:tg+tg, 1:tg+tg*2]
C33 <- C[1:tg+2*tg, 1:tg+2*tg]

beta_hat <- c(1, NA, NA)
beta_hat[2] <- sum(c(C13)*c(C23))/sum(c(C13)^2)
beta_hat[3] <- sum(c(C12)*c(C23))/sum(c(C12)^2)

coeff.C0 <- c(beta_hat[2], beta_hat[3],
              beta_hat[2]*beta_hat[3])
coeff.C1 <- beta_hat^2
C0_raw <- C1_raw <- matrix(NA, tg, tg)
for (t0 in 1:tg){
  for (t1 in 1:tg){
    C0_raw[t0, t1] <- sum(coeff.C0*c(C12[t0, t1], C13[t0, t1], 
                                     C23[t0, t1]))/sum(coeff.C0^2)
  }
}

C0_raw <- forceSymmetric(C0_raw)
C0_fit <- sFPCA_fit(C0_raw)
C0 <- C0_fit$S

for (t0 in 1:tg){
  for (t1 in 1:tg){
    C1_raw[t0, t1] <- sum(coeff.C1*c(C11[t0, t1] - beta_hat[1]^2*C0[t0, t1], 
                                     C22[t0, t1] - beta_hat[2]^2*C0[t0, t1],
                                     C33[t0, t1] - beta_hat[3]^2*C0[t0, t1]))/sum(coeff.C1^2)
  }
}

C1_fit <- sFPCA_fit(C1_raw)
C1 <- C1_fit$S

face.fit <- sFPCA_post_all(C0_fit, C1_fit, fit.mface, beta_hat, pve = 0.99)
d0_est <- face.fit$eigenvalues0[1:L0]
phi_est <- face.fit$eigenfunctions0[, 1:L0]
d1_est <- face.fit$eigenvalues1[1:L1]
psi_est <- matrix(face.fit$eigenfunctions1[, 1:L1], nrow = tg)

## Revert the sign of eigenfunctions for correct interpretation of eigenfunctions
phi_sign <- c(1, 1)
phi_index <- c(51, 1)
phi <- matrix(NA, nrow(long.train), ncol = L0)
phi_surv <- matrix(NA, N.train, ncol = L0)
for (l in 1:L0){
  phi_est[, l] <- sign_eigen(phi_est[, l], phi_index[l], phi_sign[l])
  phi.smooth <- bs.smooth(phi_est[, l], tnew, argvals.new = long.train$Time, nbasis = P)
  phi[, l] <- phi.smooth$est.value
  phi.smooth <- bs.smooth(phi_est[, l], tnew, argvals.new = surv.train$surv_time, nbasis = P)
  phi_surv[, l] <- phi.smooth$est.value
}

psi_sign <- c(1)
psi_index <- c(1)
psi <- matrix(NA, nrow(long.train), ncol = L1)
psi_surv <- matrix(NA, N.train, ncol = L1)
for (l in 1:L1){
  psi_est[, l] <- sign_eigen(psi_est[, l], psi_index[l], psi_sign[l])
  psi.smooth <- bs.smooth(psi_est[, l], tnew, argvals.new = long.train$Time, nbasis = P)
  psi[, l] <- psi.smooth$est.value
  psi.smooth <- bs.smooth(psi_est[, l], tnew, argvals.new = surv.train$surv_time, nbasis = P)
  psi_surv[, l] <- psi.smooth$est.value
}

sigma_hat <- colMeans(sqrt(face.fit$var.error.new))
bs.mat <- create_bs(time.grid = tnew, pred.time = long.train$Time, nbasis = P)
bs.mat.grid <- create_bs(tnew, tnew, nbasis = P)

## Determine missingness
Y_missing <- missing_determine(value = Y.list.o)
mean_Y <- list(mean1 = rep(mean(long.train$y1, na.rm = T), Y_missing$len.missing[1]), 
               mean2 = rep(mean(long.train$y2, na.rm = T), Y_missing$len.missing[2]),
               mean3 = rep(mean(long.train$y3, na.rm = T), Y_missing$len.missing[3]))
Y.list.new <- miss.index <- list()
for (j in 1:J){
  Y.list.new[[j]] <- Y_missing$new.value[[j]]
  miss.index[[j]] <- Y_missing$missing.index[[j]]
}

G <- 15 ## number of Gaussian quadrature points
weights <- gauss.quad(G, kind = 'legendre')$weights ## w_g^* in document
nodes <- gauss.quad(G, kind = 'legendre')$nodes  ## x_g^* in document

## Stan sampling ##
md = stan_model('source_code/M2_NM.stan')
stan_dat <- list(n = N.train, J = J, 
                 nobs = nrow(long.train), nmiss = Y_missing$len.missing, 
                 id_long = long.train$ID, 
                 L0 = L0, L1 = L1, 
                 P = P, P_surv = ncol(x), 
                 Y1 = Y.list.new[[1]], Y2 = Y.list.new[[2]], Y3 = Y.list.new[[3]],
                 miss_index1 = miss.index[[1]], miss_index2 = miss.index[[2]], 
                 miss_index3 = miss.index[[3]], 
                 time = long.train$Time, x = x, 
                 surv_time = surv.train$surv_time, status = surv.train$status,
                 b = bs.mat)

rnd <- 1/3
inits1 <- list(A1 = g.inits(mu_est[[1]]$coefficient), A2 = g.inits(mu_est[[2]]$coefficient),
               A3 = g.inits(mu_est[[3]]$coefficient), 
               sigma = g.inits.sd(sigma_hat),
               logh0 = g.inits(logh0_hat), 
               gamma_x = as.array(g.inits(gamma_x_hat)), 
               Y1_imp = mean_Y[[1]]+rnd, Y2_imp = mean_Y[[2]]+rnd, Y3_imp = mean_Y[[3]]+rnd)
inits2 <- list(A1 = g.inits(mu_est[[1]]$coefficient), A2 = g.inits(mu_est[[2]]$coefficient),
               A3 = g.inits(mu_est[[3]]$coefficient), 
               sigma = g.inits.sd(sigma_hat),
               logh0 = g.inits(logh0_hat), 
               gamma_x = as.array(g.inits(gamma_x_hat)), 
               Y1_imp = mean_Y[[1]]+rnd, Y2_imp = mean_Y[[2]]+rnd, Y3_imp = mean_Y[[3]]+rnd)

inits <- list(c1 = inits1, c2 = inits2)
pars <- c('A1', 'A2', 'A3',
          'sigma', 'logh0', 'gamma_x')
n.iters = 3000
n.warmups = 2000
fitStan <- sampling(md, data = stan_dat, iter = n.iters, warmup = n.warmups, 
                    chains = 2, thin=1, init = inits, pars = pars, seed = 2020*seed,
                    control = list(adapt_delta = 0.8, max_treedepth=10))
fname <- paste0('stan_fit/', model, '/', seed, '_posterior.RData')
save(list = 'fitStan', file = fname)

Q <- (n.iters - n.warmups)*2 ## samples
A_sample <- extract(fitStan, pars = c('A1', 'A2', 'A3'))
sigma_sample <- extract(fitStan, pars = 'sigma')[[1]]
logh0_sample <- extract(fitStan, pars = 'logh0')[[1]]
gamma_x_sample <- extract(fitStan, pars = 'gamma_x')[[1]]

mu1_hat <- mu2_hat <- mu3_hat <- matrix(0, Q, tg)
for (i in 1:Q){
  mu1_hat[i, ] <- (bs.mat.grid %*% A_sample[[1]][i, ])[, 1]
  mu2_hat[i, ] <- (bs.mat.grid %*% A_sample[[2]][i, ])[, 1]
  mu3_hat[i, ] <- (bs.mat.grid %*% A_sample[[3]][i, ])[, 1]
}

summary.parameters <- data.frame(mean = 0, sd = 0, lower_25_ci = 0, upper_975_ci = 0)

i1 <- 0
summary.parameters[(i1+1):(i1+J), ] <- get_quantile_matrix(sigma_sample)

i2 <- i1 + J + 1
summary.parameters[i2, ] <- get_quantile_vector(logh0_sample)
summary.parameters[i2+1, ] <- get_quantile_matrix(gamma_x_sample)

summary.mean <- data.frame(mean = 0)
summary.mean[1:tg, ] <- get_mean_function(mu1_hat)
summary.mean[(1+tg):(tg*2), ] <- get_mean_function(mu2_hat)
summary.mean[(1+tg*2):(tg*3), ] <- get_mean_function(mu3_hat)

fname <- paste0('result/', model, '/', seed, '_posterior.RData')
save(list = c('summary.parameters', 'summary.mean'), file = fname)

## Dynamic prediction
A <- lapply(1:J, FUN = function(x) colMeans(A_sample[[x]]))
sigma <- sapply(1:J, FUN = function(x) mean(sigma_sample[, x]))
logh0 <- mean(logh0_sample)
gamma_x <- apply(gamma_x_sample, 2, 'mean')

## compute true survival probability at time t: exp(-cumulative hazard)
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

## compute estimated survival probability at time t: exp(-cumulative hazard)
H_est <- function(x, t){
  H <- exp(logh0 + x*gamma_x)*t
  return(H)
}

starting.time <- c(0.3, 0.4, 0.5, 0.55, 0.6) 
delta.time <- seq(0.1, 0.25, by = 0.01)
AUC <- BS <- matrix(NA, length(starting.time), length(delta.time))
AUC.true <- BS.true <- matrix(NA, length(starting.time), length(delta.time))

for (Tstart in starting.time){
  
  ## filter out subjects with observed survival time greater than starting time
  surv.test2 <- surv.test[surv.test$surv_time>Tstart, ]
  long.test2 <- long.test[long.test$ID %in% surv.test2$ID, ]
  surv.test2$newID <- 1:nrow(surv.test2)
  long.test2$newID <- rep(1:nrow(surv.test2), table(long.test2$ID))
  x.test <- as.matrix(surv.test2$x)
  
  true_xi <- dat$xi[surv.test2$ID, ]
  true_zeta1 <- as.matrix(dat$zeta[[1]][surv.test2$ID, ])
  true_zeta2 <- as.matrix(dat$zeta[[2]][surv.test2$ID, ])
  true_zeta3 <- as.matrix(dat$zeta[[3]][surv.test2$ID, ])
  true_xi_m <- dat$xi_m[surv.test2$ID, ]
  
  ## filter out longitudinal observations prior or equal to starting time
  long.test.prior <- long.test2[long.test2$Time<=Tstart, ]
  long.test.posterior <- long.test2[long.test2$Time>Tstart, ]
  tmp.ID.test <- unique(long.test2$newID)
  bs.mat.test <- create_bs(time.grid = tnew, pred.time = long.test.prior$Time, nbasis = P)
  
  ## Determine missingness
  Y.list.test <- list(y1 = long.test.prior$y1, y2 = long.test.prior$y2, y3 = long.test.prior$y3)
  Y_missing.test <- missing_determine(value = Y.list.test)
  
  ## calculate conditional survival probabilities at T+\delta_T
  for (Tdelta in delta.time){
    tmp.surv.true <- rep(NA, length(tmp.ID.test))
    tmp.surv.predict <- rep(NA, length(tmp.ID.test))
    
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
      
      tmp.surv.true[i] <- exp(-Hi_true_t1 + Hi_true_t0)
      
      ## estimated survival probability ##
      Hi_est_t0 <- Hi_est_t1 <- rep(NA, Q)
      for (q in 1:Q){
        Hi_est_t0[q] <- H_est(x.test[i, ], t = Tstart)
        Hi_est_t1[q] <- H_est(x.test[i, ], t = Tstart + Tdelta)
      }
      
      cond_surv_prob <- exp(-Hi_est_t1 + Hi_est_t0)
      tmp.surv.predict[i] <- mean(cond_surv_prob)
    }
    
    ROC.est <- tdROC(X = 1 - tmp.surv.predict, Y = surv.test2$surv_time,
                     delta = surv.test2$status, tau = Tstart + Tdelta,
                     span = 0.1, alpha = 0.05,
                     n.grid = 1000, cut.off = 0.5)
    AUC[which(Tstart==starting.time), which(Tdelta==delta.time)] <- ROC.est$AUC$value
    
    ROC.true <- tdROC(X = 1 - tmp.surv.true, Y = surv.test2$surv_time, 
                      delta = surv.test2$status, tau = Tstart + Tdelta,
                      span = 0.1, alpha = 0.05,
                      n.grid = 1000, cut.off = 0.5)
    AUC.true[which(Tstart==starting.time), which(Tdelta==delta.time)] <- ROC.true$AUC$value
    
    surv.obj <- Surv(surv.test2$surv_time, surv.test2$status)
    BS[which(Tstart==starting.time), which(Tdelta==delta.time)] <- sbrier(surv.obj, tmp.surv.predict, btime = Tstart + Tdelta)
    BS.true[which(Tstart==starting.time), which(Tdelta==delta.time)] <- sbrier(surv.obj, tmp.surv.true, btime = Tstart + Tdelta)
  }
}

fname <- paste0('result/', model, '/', seed, '_dp.RData')
save(list = c('AUC', 'BS', 'AUC.true', 'BS.true'), file = fname)

