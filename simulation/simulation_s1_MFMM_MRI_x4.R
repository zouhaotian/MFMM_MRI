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

model = 's1_MFMM_MRI_x4'

source('source_code/functions.R')
source('source_code/data_gen1.R')
source('source_code/sFPCA.R')
source('source_code/search.R')
load('RData/data_genx4.RData')

# seed <- 1
N.train <- tp$N.train; N.test <- tp$N.test; N <- N.train + N.test
ID <- 1:N; ID.train <- 1:N.train; ID.test <- (N.train+1):N

dat <- sim(seed*2019, tp, N)
long <- dat$long; long.train <- long[long$ID %in% ID.train, ]; long.test <- long[long$ID %in% ID.test, ]
surv <- dat$surv; surv.train <- surv[surv$ID %in% ID.train, ]; surv.test <- surv[surv$ID %in% ID.test, ]
MRI.dat = dat$MRI.dat; M.train <- MRI.dat[ID.train, ]; M.test <- MRI.dat[ID.test, ]
xi_m_true = dat$xi_m; xi_m_true.train = xi_m_true[ID.train, ]; xi_m_true.test = xi_m_true[ID.test, ]
xi_true = dat$xi; xi_true.train = xi_true[ID.train, ]; xi_true.test = xi_true[ID.test, ]
C0_true = dat$C0; C1_true = dat$C1; Ct_true = dat$Ct

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

#####

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

## weights; nodes; time vector of length N denoting interval (0, time[i]) ##
create_gq <- function(w, x, time){
  const <- rep(NA, length(time)) ## const_{i} = (b_{i} - a_{i})/2
  a <- b <- rep(NA, length(time))
  
  phi1_interval <- phi2_interval <- array(NA, dim = c(length(time), G))
  psi1_interval <- array(NA, dim = c(length(time), G))
  for (i in 1:length(time)){
    a[i] <- 0
    b[i] <- time[i]
    const[i] <- (b[i] - a[i])/2
    x_star <- (b[i] - a[i])/2*x + (b[i] + a[i])/2
    
    ## Evaluate interval of phi and psi on x_star
    phi.smooth <- bs.smooth(phi_est[, 1], tnew, argvals.new = x_star, nbasis = P)
    phi1_interval[i, ] <- phi.smooth$est.value
    phi.smooth <- bs.smooth(phi_est[, 2], tnew, argvals.new = x_star, nbasis = P)
    phi2_interval[i, ] <- phi.smooth$est.value
    
    psi.smooth <- bs.smooth(psi_est[, 1], tnew, argvals.new = x_star, nbasis = P)
    psi1_interval[i, ] <- psi.smooth$est.value
  }
  l <- list(const = const, 
            phi1_interval = phi1_interval, 
            phi2_interval = phi2_interval,
            psi1_interval = psi1_interval)
  return(l)
}

## Create Gaussian quadrature points on grid (a_{i}, b_{i}) ##
G <- 15 ## number of Gaussian quadrature points
weights <- gauss.quad(G, kind = 'legendre')$weights ## w_g^* in document
nodes <- gauss.quad(G, kind = 'legendre')$nodes  ## x_g^* in document
l <- create_gq(weights, nodes, surv.train$surv_time)
const <- l$const
phi1_interval <- l$phi1_interval
phi2_interval <- l$phi2_interval
psi1_interval <- l$psi1_interval

## Stan sampling ##
md = stan_model('source_code/M2_prior2.stan')
stan_dat <- list(n = N.train, J = J, 
                 nobs = nrow(long.train), nmiss = Y_missing$len.missing, 
                 id_long = long.train$ID, 
                 L0 = L0, L1 = L1, 
                 P = P, P_surv = ncol(x), 
                 Y1 = Y.list.new[[1]], Y2 = Y.list.new[[2]], Y3 = Y.list.new[[3]],
                 miss_index1 = miss.index[[1]], miss_index2 = miss.index[[2]], 
                 miss_index3 = miss.index[[3]], 
                 time = long.train$Time, x = x, 
                 b = bs.mat, phi = phi, psi = psi)

rnd <- 1/3
inits1 <- list(A1 = g.inits(mu_est[[1]]$coefficient), A2 = g.inits(mu_est[[2]]$coefficient),
               A3 = g.inits(mu_est[[3]]$coefficient), 
               sqrt_d0 = as.array(g.inits.sd(sqrt(d0_est))), 
               sqrt_d1 = as.array(g.inits.sd(sqrt(d1_est))),
               beta = g.inits(beta_hat[2:J]), 
               sigma = g.inits.sd(sigma_hat),
               Y1_imp = mean_Y[[1]]+rnd, Y2_imp = mean_Y[[2]]+rnd, Y3_imp = mean_Y[[3]]+rnd)
inits2 <- list(A1 = g.inits(mu_est[[1]]$coefficient), A2 = g.inits(mu_est[[2]]$coefficient),
               A3 = g.inits(mu_est[[3]]$coefficient), 
               sqrt_d0 = as.array(g.inits.sd(sqrt(d0_est))), 
               sqrt_d1 = as.array(g.inits.sd(sqrt(d1_est))),
               beta = g.inits(beta_hat[2:J]), 
               sigma = g.inits.sd(sigma_hat),
               Y1_imp = mean_Y[[1]]+rnd, Y2_imp = mean_Y[[2]]+rnd, Y3_imp = mean_Y[[3]]+rnd)

inits <- list(c1 = inits1, c2 = inits2)
pars <- c('A1', 'A2', 'A3', 'd0', 'd1', 
          'beta', 'sigma', 'xi')
n.iters = 3000
n.warmups = 2000
fitStan <- sampling(md, data = stan_dat, iter = n.iters, warmup = n.warmups, 
                    chains = 2, thin=1, init = inits, pars = pars, seed = 2020*seed,
                    control = list(adapt_delta = 0.8, max_treedepth=10))
fname <- paste0('stan_fit/', model, '/', seed, '_prior2.RData')
save(list = 'fitStan', file = fname)
# load(fname)

A_sample <- extract(fitStan, pars = c('A1', 'A2', 'A3'))
A = lapply(1:J, FUN = function(x) colMeans(A_sample[[x]]))

beta_sample = extract(fitStan, pars = 'beta')[[1]]
beta = c(1, colMeans(beta_sample))

d0_sample = extract(fitStan, pars = 'd0')[[1]]
d0 = colMeans(d0_sample)

d1_sample = extract(fitStan, pars = 'd1')[[1]]
d1 = colMeans(d1_sample)

sigma_sample = extract(fitStan, pars = 'sigma')[[1]]
sigma = colMeans(sigma_sample)

xi_sample <- extract(fitStan, pars = 'xi')[[1]] 
summary.xi = t(colMeans(xi_sample)) ## is a (N*L0) matrix
summary.xi = cbind(summary.xi, xi_true.train)

## Filter out longitudinal data with more observations available
nobs = as.vector(table(long.train$ID))
dat2 = data.frame(ID = unique(long.train$ID), nobs = nobs)
dat3 = dat2[which(dat2$nobs>=7), ]   ## We can try different numbers of observation
tmp.ID = dat3$ID

## 
Pm = round(V/30)
mu_m = colMeans(M.train)
mu_m_smooth = bs.smooth(mu_m, obsgrid, obsgrid, nbasis = Pm)
mu_m_expand = matrix(rep(mu_m_smooth$est.value, N.train), nrow = N.train, byrow = T)
M.demean = M.train - mu_m_expand

phi_m.sign = c(1, 1)
phi_m.index = c(V, round(V/3))
psi_m.sign = rep(1, Lm)
psi_m.index = rep(1, Lm)

sign_beta_m = -1

## Estimate eigenfunction of m_i(v) ##
ID1 = tmp.ID
ID2 = tmp.ID
xi1 = summary.xi[ID1, 1:L0]; xi2 = summary.xi[ID2, 1:L0]
#xi1 = xi_true.train[ID1, ]; xi2 = xi_true.train[ID2, ]
M1 = M.demean[ID1, ]; M2 = M.demean[ID2, ]

res = search(M1, M2, xi1, xi2, L0, Lm, V, 
             phi_m.sign, phi_m.index, 
             psi_m.sign, psi_m.index,
             max.iters = 20, eps = 1e-5, 
             beta_m = NULL, Phi_m_est = NULL,
             d0 = d0, sign_beta_m)
beta_m <- res$beta_m
dm <- res$dm
Phi_m_est <- res$Phi_m_est
Psi_m_est <- res$Psi_m_est
sigma_m = res$sigma_m

## Create m_{ij} matrix ##
m_mat = matrix(NA, nrow = N.train, ncol = Lm)
f_l = matrix(NA, nrow = L0, ncol = Lm)

for (j in 1:Lm){
  for (i in 1:N.train){
    m_mat[i, j] <- c_int(M.demean[i, ], Psi_m_est[, j], width = 1/(V-1))
  }
  for (l in 1:L0){
    f_l[l, j] <- c_int(Phi_m_est[, l], Psi_m_est[, j], width = 1/(V-1))
  }
}

## 
md = stan_model('source_code/M2.stan')
stan_dat <- list(n = N.train, J = J, 
                 nobs = nrow(long.train), nmiss = Y_missing$len.missing, 
                 id_long = long.train$ID, 
                 L0 = L0, L1 = L1, Lm = Lm, 
                 P = P, P_surv = P_surv, 
                 Y1 = Y.list.new[[1]], Y2 = Y.list.new[[2]], Y3 = Y.list.new[[3]],
                 miss_index1 = miss.index[[1]], miss_index2 = miss.index[[2]], 
                 miss_index3 = miss.index[[3]], 
                 time = long.train$Time,
                 x = x, surv_time = surv.train$surv_time, status = surv.train$status,
                 b = bs.mat, phi = phi, psi = psi, 
                 phi_surv = phi_surv, psi_surv = psi_surv,
                 m_mat = m_mat, f_l = f_l, beta_m = beta_m, 
                 G = G, w = weights, constant = const, 
                 phi1_interval = phi1_interval, phi2_interval = phi2_interval,
                 psi1_interval = psi1_interval)
rnd <- 1/3
inits1 <- list(A1 = g.inits(A[[1]]), A2 = g.inits(A[[2]]), A3 = g.inits(A[[3]]), 
               sqrt_d0 = as.array(g.inits.sd(sqrt(d0))), 
               sqrt_d1 = as.array(g.inits.sd(sqrt(d1))),
               xi = t(summary.xi[, 1:L0]), 
               beta = g.inits(beta[2:J]), 
               sigma = g.inits.sd(sigma),
               logh0 = g.inits(logh0_hat), 
               gamma_x = as.array(g.inits(gamma_x_hat)), 
               gamma0 = rnorm(1, 0, 1/3),
               gamma1 = as.array(rnorm(J, 0, 1/3)), 
               gamma_m = as.array(rnorm(Lm, 0, 1/3)), 
               Y1_imp = mean_Y[[1]]+rnd, Y2_imp = mean_Y[[2]]+rnd, Y3_imp = mean_Y[[3]]+rnd)
inits2 <- list(A1 = g.inits(A[[1]]), A2 = g.inits(A[[2]]), A3 = g.inits(A[[3]]), 
               sqrt_d0 = as.array(g.inits.sd(sqrt(d0))), 
               sqrt_d1 = as.array(g.inits.sd(sqrt(d1))),
               xi = t(summary.xi[, 1:L0]), 
               beta = g.inits(beta[2:J]), 
               sigma = g.inits.sd(sigma),
               logh0 = g.inits(logh0_hat), 
               gamma_x = as.array(g.inits(gamma_x_hat)), 
               gamma0 = rnorm(1, 0, 1/3),
               gamma1 = as.array(rnorm(J, 0, 1/3)), 
               gamma_m = as.array(rnorm(Lm, 0, 1/3)), 
               Y1_imp = mean_Y[[1]]-rnd, Y2_imp = mean_Y[[2]]-rnd, Y3_imp = mean_Y[[3]]-rnd)
inits <- list(c1 = inits1, c2 = inits2)
pars <- c('A1', 'A2', 'A3', 'd0', 'd1', 'dm',
          'beta',
          'sigma',
          'logh0', 'gamma_x', 'gamma0', 
          'gamma1', 'gamma_m', 
          'xi', 'xi_m')
n.iters = 5000
n.warmups = 4000
fitStan <- sampling(md, data = stan_dat, iter = n.iters, warmup = n.warmups, 
                    chains = 2, thin=1, init = inits, pars = pars, seed = 2020*seed,
                    control = list(adapt_delta = 0.8, max_treedepth=10))
fname <- paste0('stan_fit/', model, '/', seed, '_posterior.RData')
save(list = 'fitStan', file = fname)
# load(fname)

Q <- (n.iters - n.warmups)*2 ## samples
A_sample <- extract(fitStan, pars = c('A1', 'A2', 'A3'))
d0_sample <- extract(fitStan, pars = 'd0')[[1]]
d1_sample <- extract(fitStan, pars = 'd1')[[1]]
dm_sample <- extract(fitStan, pars = 'dm')[[1]]
beta_sample <- extract(fitStan, pars = 'beta')[[1]]
sigma_sample <- extract(fitStan, pars = 'sigma')[[1]]
logh0_sample <- extract(fitStan, pars = 'logh0')[[1]]
gamma_x_sample <- extract(fitStan, pars = 'gamma_x')[[1]]
gamma0_sample <- extract(fitStan, pars = 'gamma0')[[1]]
gamma1_sample <- extract(fitStan, pars = 'gamma1')[[1]]
gamma_m_sample <- extract(fitStan, pars = 'gamma_m')[[1]]

mu1_hat <- mu2_hat <- mu3_hat <- matrix(0, Q, tg)
for (i in 1:Q){
  mu1_hat[i, ] <- (bs.mat.grid %*% A_sample[[1]][i, ])[, 1]
  mu2_hat[i, ] <- (bs.mat.grid %*% A_sample[[2]][i, ])[, 1]
  mu3_hat[i, ] <- (bs.mat.grid %*% A_sample[[3]][i, ])[, 1]
}

summary.parameters <- data.frame(mean = 0, sd = 0, lower_25_ci = 0, upper_975_ci = 0)
summary.parameters[1:L0, ] <- get_quantile_matrix(d0_sample)
summary.parameters[(L0+1):(L0+L1), ] <- get_quantile_matrix(d1_sample)
summary.parameters[(L0+L1+1):(L0+L1+Lm), ] <- get_quantile_matrix(dm_sample)

i1 <- L0 + L1 + Lm
summary.parameters[(i1+1):(i1+J-1), ] <- get_quantile_matrix(beta_sample)
summary.parameters[(i1+J):(i1+J), c(1, 3, 4)] <- beta_m
summary.parameters[(i1+J+1):(i1+J*2), ] <- get_quantile_matrix(sigma_sample)
summary.parameters[(i1+J*2+1):(i1+J*2+1), ] <- sigma_m

i2 <- i1 + J*2 + 2
summary.parameters[i2, ] <- get_quantile_vector(logh0_sample)
summary.parameters[i2+1, ] <- get_quantile_matrix(gamma_x_sample)
summary.parameters[i2+2, ] <- get_quantile_vector(gamma0_sample)

i3 <- i2 + 2
summary.parameters[(i3+1):(i3+J), ] <- get_quantile_matrix(gamma1_sample)

i4 = i3 + J
summary.parameters[(i4+1):(i4+Lm), ] <- get_quantile_matrix(gamma_m_sample)

summary.mean <- data.frame(mean = 0)
summary.mean[1:tg, ] <- get_mean_function(mu1_hat)
summary.mean[(1+tg):(tg*2), ] <- get_mean_function(mu2_hat)
summary.mean[(1+tg*2):(tg*3), ] <- get_mean_function(mu3_hat)
summary.mean[(1+tg*3):(tg*3+V), ] <- mu_m_smooth$est.value

xi_sample <- extract(fitStan, pars = 'xi')[[1]] 
summary.xi = t(colMeans(xi_sample)) ## is a (N*L0) matrix
summary.xi = cbind(summary.xi, xi_true.train)

xi_m_sample <- extract(fitStan, pars = 'xi_m')[[1]] 
summary.xi_m = colMeans(xi_m_sample) ## is a (N*Lm) matrix
summary.xi_m = cbind(summary.xi_m, xi_m_true.train)

fname <- paste0('result/', model, '/', seed, '_posterior.RData')
save(list = c('summary.parameters', 'summary.mean', 'summary.xi', 'summary.xi_m', 
              'phi_est', 'psi_est', 'Phi_m_est', 'Psi_m_est'), file = fname)

## Dynamic prediction
A <- lapply(1:J, FUN = function(x) colMeans(A_sample[[x]]))
d0 <- apply(d0_sample, 2, 'mean')
d1 <- apply(d1_sample, 2, 'mean')
beta <- apply(beta_sample, 2, 'mean')
sigma <- sapply(1:J, FUN = function(x) mean(sigma_sample[, x]))
logh0 <- mean(logh0_sample)
gamma_x <- apply(gamma_x_sample, 2, 'mean')
gamma0 <- mean(gamma0_sample)
gamma1 <- apply(gamma1_sample, 2, 'mean')
gamma_m <- apply(gamma_m_sample, 2, 'mean')

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
H_est <- function(phi1_interval_t, phi2_interval_t, psi1_interval_t,
                  x, weights, xi, zeta1, zeta2, zeta3, xi_m, t){
  U = xi[1]*phi1_interval_t + xi[2]*phi2_interval_t
  W1 = zeta1*psi1_interval_t
  W2 = zeta2*psi1_interval_t
  W3 = zeta3*psi1_interval_t
  f = exp(gamma0*U + gamma1[1]*W1 + gamma1[2]*W2 + gamma1[3]*W3)
  H <- exp(logh0 + x*gamma_x + sum(gamma_m*xi_m))*t/2*sum(f*weights)
  return(H)
}

md_dp <- stan_model('source_code/M2_dp.stan')
pars_dp <- c('xi', 'zeta_1', 'zeta_2', 'zeta_3', 'xi_m')
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
  M.test2 <- M.test[surv.test2$ID - N.train, ]
  mu_m_expand = matrix(rep(mu_m_smooth$est.value, nrow(surv.test2)), nrow = nrow(surv.test2), byrow = T)
  M.demean = M.test2 - mu_m_expand
  
  true_xi <- dat$xi[surv.test2$ID, ]
  true_zeta1 <- as.matrix(dat$zeta[[1]][surv.test2$ID, ])
  true_zeta2 <- as.matrix(dat$zeta[[2]][surv.test2$ID, ])
  true_zeta3 <- as.matrix(dat$zeta[[3]][surv.test2$ID, ])
  true_xi_m <- dat$xi_m[surv.test2$ID, ]
  
  ## Create m_{ij} matrix ##
  m_mat = matrix(NA, nrow = nrow(surv.test2), ncol = Lm)
  f_l = matrix(NA, nrow = L0, ncol = Lm)
  
  for (j in 1:Lm){
    for (i in 1:nrow(surv.test2)){
      m_mat[i, j] <- c_int(M.demean[i, ], Psi_m_est[, j], width = 1/(V-1))
    }
    for (l in 1:L0){
      f_l[l, j] <- c_int(Phi_m_est[, l], Psi_m_est[, j], width = 1/(V-1))
    }
  }
  
  ## filter out longitudinal observations prior or equal to starting time
  long.test.prior <- long.test2[long.test2$Time<=Tstart, ]
  long.test.posterior <- long.test2[long.test2$Time>Tstart, ]
  tmp.ID.test <- unique(long.test2$newID)
  bs.mat.test <- create_bs(time.grid = tnew, pred.time = long.test.prior$Time, nbasis = P)
  
  phi_test <- matrix(NA, nrow(long.test.prior), ncol = L0)
  for (l in 1:L0){
    phi.smooth <- bs.smooth(phi_est[, l], tnew, argvals.new = long.test.prior$Time, nbasis = P)
    phi_test[, l] <- phi.smooth$est.value
  }
  
  psi_test <- matrix(NA, nrow(long.test.prior), ncol = L1)
  for (l in 1:L1){
    psi.smooth <- bs.smooth(psi_est[, l], tnew, argvals.new = long.test.prior$Time, nbasis = P)
    psi_test[, l] <- psi.smooth$est.value
  }
  
  ## Determine missingness
  Y.list.test <- list(y1 = long.test.prior$y1, y2 = long.test.prior$y2, y3 = long.test.prior$y3)
  Y_missing.test <- missing_determine(value = Y.list.test)
  
  ## Create Gaussian quadrature points on grid (a_{i}, b_{i}) ##
  l_T0 <- create_gq(weights, nodes, time = Tstart)
  const_T0 <- l_T0$const
  phi1_interval_T0 <- l_T0$phi1_interval
  phi2_interval_T0 <- l_T0$phi2_interval
  psi1_interval_T0 <- l_T0$psi1_interval
  
  ## perform NUTS for posterior sampling of xi and zeta 
  stan_dat_dp <- list(n = length(tmp.ID.test), J = J, 
                      nobs = nrow(long.test.prior),
                      id_long = long.test.prior$newID,
                      L0 = L0, L1 = L1, Lm = Lm, 
                      P = P, P_surv = ncol(x.test), 
                      Y1 = Y_missing.test$new.value[[1]], 
                      Y2 = Y_missing.test$new.value[[2]], 
                      Y3 = Y_missing.test$new.value[[3]],
                      time = long.test.prior$Time, x = x.test,
                      surv_time = rep(Tstart, length(tmp.ID.test)), 
                      b = bs.mat.test, phi = phi_test, psi = psi_test,
                      m_mat = m_mat, f_l = f_l, beta_m = beta_m, 
                      G = G, w = weights, constant = const_T0, 
                      phi1_interval = phi1_interval_T0, phi2_interval = phi2_interval_T0,
                      psi1_interval = psi1_interval_T0,
                      A1 = A[[1]], A2 = A[[2]], A3 = A[[3]], 
                      sqrt_d0 = as.array(sqrt(d0)), sqrt_d1 = as.array(sqrt(d1)), 
                      beta = beta, sigma = sigma,
                      logh0 = logh0, gamma_x = as.array(gamma_x), 
                      gamma0 = gamma0, 
                      gamma1 = gamma1, 
                      gamma_m = as.array(gamma_m))
  inits1_dp <- inits2_dp <- list(xi = matrix(0, L0, length(tmp.ID.test)), 
                                 zeta1 = matrix(0, L1, length(tmp.ID.test)), 
                                 zeta2 = matrix(0, L1, length(tmp.ID.test)), 
                                 zeta3 = matrix(0, L1, length(tmp.ID.test)))
  inits_dp <- list(c1 = inits1_dp, c2 = inits2_dp)
  fitStan_dp <- sampling(md_dp, data = stan_dat_dp, iter = 2000, warmup = 1000, 
                         chains = 2, thin=1, init = inits_dp, pars = pars_dp, seed = 123,
                         control = list(adapt_delta = 0.8, max_treedepth=10))
  xi_sample <- extract(fitStan_dp, pars = 'xi')$xi ## extract xi_{il}
  zeta1_sample <- extract(fitStan_dp, pars = 'zeta_1')$zeta_1 ## extract zeta1_{il}
  zeta2_sample <- extract(fitStan_dp, pars = 'zeta_2')$zeta_2 ## extract zeta2_{il}
  zeta3_sample <- extract(fitStan_dp, pars = 'zeta_3')$zeta_3
  xi_m_sample <- extract(fitStan_dp, pars = 'xi_m')$xi_m ## extract xi_{il}
  
  ## calculate conditional survival probabilities at T+\delta_T
  for (Tdelta in delta.time){
    l_T1 <- create_gq(weights, nodes, time = Tstart+Tdelta)
    const_T1 <- l_T1$const
    phi1_interval_T1 <- l_T1$phi1_interval[1, ]
    phi2_interval_T1 <- l_T1$phi2_interval[1, ]
    psi1_interval_T1 <- l_T1$psi1_interval[1, ]
    
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
      xi_i <- xi_sample[, , i] 
      zeta1_i <- zeta1_sample[, , i]
      zeta2_i <- zeta2_sample[, , i]
      zeta3_i <- zeta3_sample[, , i]
      xi_m_i <- xi_m_sample[, i, ]
      Hi_est_t0 <- Hi_est_t1 <- rep(NA, Q)
      for (q in 1:Q){
        Hi_est_t0[q] <- H_est(phi1_interval_T0, phi2_interval_T0,
                              psi1_interval_T0, x.test[i, ], 
                              weights = weights, 
                              xi_i[q, ], zeta1_i[q], zeta2_i[q], zeta3_i[q], 
                              xi_m_i[q, ], t = Tstart)
        Hi_est_t1[q] <- H_est(phi1_interval_T1, phi2_interval_T1,
                              psi1_interval_T1, x.test[i, ], 
                              weights = weights, 
                              xi_i[q, ], zeta1_i[q], zeta2_i[q], zeta3_i[q], 
                              xi_m_i[q, ], t = Tstart + Tdelta)
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

