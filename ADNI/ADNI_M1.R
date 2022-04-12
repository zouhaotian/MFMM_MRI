#setwd("/MMFPCA_MRI/ADNI/")
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
library(loo)
rstan_options(auto_write = TRUE)

source('source_code/functions.R')
source('source_code/sFPCA.R')
source('source_code/search.R')

long <- read.csv('dataset/ADNI_long_2.csv')
surv <- read.csv('dataset/ADNI_surv.csv')
MRI.dat <- fread("dataset/hippo.mri.dat.csv") %>% as.matrix()

# M.train = matrix(NA, nrow = nrow(MRI.dat), ncol = ncol(MRI.dat))
# for (i in 1:ncol(MRI.dat)){
#  M.train[, i] = scale(MRI.dat[, i], center = F)
# }
M.train = MRI.dat

N <- nrow(surv)
J <- 5 ## number of longitudinal outcomes
V = ncol(MRI.dat)
obsgrid = (1:V)/V
P <- 9 ## number of basis functions for B-spline
L0 <- 2 ## npc for U_i(t)
L1 <- 1 ## npc for W_ij(t)

## Refine the time grid ##
tnew <- (0:101)/10
tg <- length(tnew)

fit.cox <- coxph(Surv(surv_time, status) ~ Age + Gender + Education + APOE4,
                 data = surv, method = 'breslow')
gamma_x_hat <- unname(fit.cox$coefficients)
baseline.cumulative.hazard <- basehaz(fit.cox)
x <- as.matrix(surv[, 4:7])

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

## Use mfaces package to calculate covariance matrices ##
dat.mface <- list('y1' = data.frame('subj' = long$ID, 'argvals' = long$Time, 'y' = Y.list.i[[1]]),
                  'y2' = data.frame('subj' = long$ID, 'argvals' = long$Time, 'y' = Y.list.i[[2]]),
                  'y3' = data.frame('subj' = long$ID, 'argvals' = long$Time, 'y' = Y.list.i[[3]]),
                  'y4' = data.frame('subj' = long$ID, 'argvals' = long$Time, 'y' = Y.list.i[[4]]),
                  'y5' = data.frame('subj' = long$ID, 'argvals' = long$Time, 'y' = Y.list.i[[5]]))
fit.mface <- mface.sparse(dat.mface, argvals.new = tnew, knots = 6, newdata = dat.mface)
# save(list = 'fit.mface', file = 'RData/fit_mface.RData')

## Smooth the estimated mean function ##
mu_est <- list()
for (i in 1:J){
  tmp.mu <- fit.mface$fit[[i]]$mu.new
  l <- bs.smooth(tmp.mu, tnew, tnew, nbasis = P)
  mu_est[[i]] <- list(value = tmp.mu, 
                      argvals = tnew, 
                      coefficient = l$coef.est)
}

C <- as.matrix(fit.mface$Chat.new)
Cm <- list(
  C11 = C[1:tg, 1:tg],
  C12 = C[1:tg, 1:tg+tg],
  C13 = C[1:tg, 1:tg+tg*2],
  C14 = C[1:tg, 1:tg+tg*3],
  C15 = C[1:tg, 1:tg+tg*4],
  C22 = C[1:tg+tg, 1:tg+tg],
  C23 = C[1:tg+tg, 1:tg+tg*2],
  C24 = C[1:tg+tg, 1:tg+tg*3],
  C25 = C[1:tg+tg, 1:tg+tg*4],
  C33 = C[1:tg+tg*2, 1:tg+tg*2],
  C34 = C[1:tg+tg*2, 1:tg+tg*3],
  C35 = C[1:tg+tg*2, 1:tg+tg*4],
  C44 = C[1:tg+tg*3, 1:tg+tg*3],
  C45 = C[1:tg+tg*3, 1:tg+tg*4],
  C55 = C[1:tg+tg*4, 1:tg+tg*4]
)

beta_hat <- c(1, rep(NA, J-1))
beta_hat[2] <- sum(c(Cm$C13, Cm$C14, Cm$C15)*c(Cm$C23, Cm$C24, Cm$C25))/sum(c(Cm$C13, Cm$C14, Cm$C15)^2)
beta_hat[3] <- sum(c(Cm$C12, Cm$C14, Cm$C15)*c(Cm$C23, Cm$C34, Cm$C35))/sum(c(Cm$C12, Cm$C14, Cm$C15)^2)
beta_hat[4] <- sum(c(Cm$C12, Cm$C13, Cm$C15)*c(Cm$C24, Cm$C34, Cm$C45))/sum(c(Cm$C12, Cm$C13, Cm$C15)^2)
beta_hat[5] <- sum(c(Cm$C12, Cm$C13, Cm$C14)*c(Cm$C25, Cm$C35, Cm$C45))/sum(c(Cm$C12, Cm$C13, Cm$C14)^2)

coeff.C0 <- c(beta_hat[2], beta_hat[3], beta_hat[4], beta_hat[5],
              beta_hat[2]*beta_hat[3], beta_hat[2]*beta_hat[4], beta_hat[2]*beta_hat[5],
              beta_hat[3]*beta_hat[4], beta_hat[3]*beta_hat[5],
              beta_hat[4]*beta_hat[5])
coeff.C1 <- beta_hat^2
C0_raw <- C1_raw <- matrix(NA, tg, tg)
for (t0 in 1:tg){
  for (t1 in 1:tg){
    C0_raw[t0, t1] <- sum(coeff.C0*c(Cm$C12[t0, t1], Cm$C13[t0, t1], Cm$C14[t0, t1], Cm$C15[t0, t1],
                                     Cm$C23[t0, t1], Cm$C24[t0, t1], Cm$C25[t0, t1],
                                     Cm$C34[t0, t1], Cm$C35[t0, t1],
                                     Cm$C45[t0, t1]))/sum(coeff.C0^2)
  }
}

C0_raw <- forceSymmetric(C0_raw)
C0_fit <- sFPCA_fit(C0_raw)
C0 <- C0_fit$S

for (t0 in 1:tg){
  for (t1 in 1:tg){
    C1_raw[t0, t1] <- sum(coeff.C1*c(Cm$C11[t0, t1] - beta_hat[1]^2*C0[t0, t1], 
                                     Cm$C22[t0, t1] - beta_hat[2]^2*C0[t0, t1],
                                     Cm$C33[t0, t1] - beta_hat[3]^2*C0[t0, t1],
                                     Cm$C44[t0, t1] - beta_hat[4]^2*C0[t0, t1],
                                     Cm$C55[t0, t1] - beta_hat[5]^2*C0[t0, t1]))/sum(coeff.C1^2)
  }
}

C1_fit <- sFPCA_fit(C1_raw)
C1 <- C1_fit$S

face.fit <- sFPCA_post_all(C0_fit, C1_fit, mfaces_fit = fit.mface, Beta = beta_hat, pve = 0.99)
d0_est <- face.fit$eigenvalues0[1:L0]
phi_est <- matrix(face.fit$eigenfunctions0[, 1:L0], nrow = tg)
d1_est <- face.fit$eigenvalues1[1:L1]
psi_est <- matrix(face.fit$eigenfunctions1[, 1:L1], nrow = tg)

## Revert the sign of eigenfunctions for correct interpretation of eigenfunctions
phi_sign <- c(1, 1)
phi_index <- c(1, 1)
phi <- matrix(NA, nrow(long), ncol = L0)
for (l in 1:L0){
  phi_est[, l] <- sign_eigen(phi_est[, l], phi_index[l], phi_sign[l])
  phi.smooth <- bs.smooth(phi_est[, l], tnew, argvals.new = long$Time, nbasis = P)
  phi[, l] <- phi.smooth$est.value
}

psi_sign <- c(1)
psi_index <- c(1)
psi <- matrix(NA, nrow(long), ncol = L1)
for (l in 1:L1){
  psi_est[, l] <- sign_eigen(psi_est[, l], psi_index[l], psi_sign[l])
  psi.smooth <- bs.smooth(psi_est[, l], tnew, argvals.new = long$Time, nbasis = P)
  psi[, l] <- psi.smooth$est.value
}

#save(list = c('tnew', 'phi_est', 'psi_est'), file ='RData/est_efunc.RData')
sigma_hat <- colMeans(sqrt(face.fit$var.error.new))
bs.mat <- create_bs(time.grid = tnew, pred.time = long$Time, nbasis = P)

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

## Stan fit ##
md = stan_model('source_code/M1_prior2.stan')
stan_dat <- list(n = N, J = J, 
                 nobs = nrow(long), nmiss = Y_missing$len.missing, 
                 id_long = long$ID, 
                 L0 = L0, L1 = L1, P = P,
                 Y1 = Y.list.new[[1]], Y2 = Y.list.new[[2]], Y3 = Y.list.new[[3]],
                 Y4 = Y.list.new[[4]], Y5 = Y.list.new[[5]],
                 miss_index1 = miss.index[[1]], miss_index2 = miss.index[[2]], 
                 miss_index3 = miss.index[[3]], miss_index4 = miss.index[[4]], 
                 miss_index5 = miss.index[[5]],
                 time = long$Time,
                 b = bs.mat, phi = phi, psi = psi)

rnd <- 1/3
inits1 <- list(A1 = g.inits(mu_est[[1]]$coefficient), A2 = g.inits(mu_est[[2]]$coefficient),
               A3 = g.inits(mu_est[[3]]$coefficient), A4 = g.inits(mu_est[[4]]$coefficient), 
               A5 = g.inits(mu_est[[5]]$coefficient), 
               sqrt_d0 = as.array(g.inits.sd(sqrt(d0_est))), 
               sqrt_d1 = as.array(g.inits.sd(sqrt(d1_est))),
               beta = g.inits(beta_hat[2:J]), 
               sigma = g.inits.sd(sigma_hat),
               Y1_imp = mean_Y[[1]]+rnd, Y2_imp = mean_Y[[2]]+rnd, Y3_imp = mean_Y[[3]]+rnd, 
               Y4_imp = mean_Y[[4]]+rnd, Y5_imp = mean_Y[[5]]+rnd)
inits2 <- list(A1 = g.inits(mu_est[[1]]$coefficient), A2 = g.inits(mu_est[[2]]$coefficient),
               A3 = g.inits(mu_est[[3]]$coefficient), A4 = g.inits(mu_est[[4]]$coefficient), 
               A5 = g.inits(mu_est[[5]]$coefficient), 
               sqrt_d0 = as.array(g.inits.sd(sqrt(d0_est))), 
               sqrt_d1 = as.array(g.inits.sd(sqrt(d1_est))),
               beta = g.inits(beta_hat[2:J]), 
               sigma = g.inits.sd(sigma_hat),
               Y1_imp = mean_Y[[1]]-rnd, Y2_imp = mean_Y[[2]]-rnd, Y3_imp = mean_Y[[3]]-rnd, 
               Y4_imp = mean_Y[[4]]-rnd, Y5_imp = mean_Y[[5]]-rnd)

inits <- list(c1 = inits1, c2 = inits2)
pars <- c('A1', 'A2', 'A3', 'A4', 'A5', 
          'd0', 'd1', 'beta', 
          'sigma', 
          'xi', 'zeta_1', 'zeta_2', 'zeta_3', 'zeta_4', 'zeta_5',
          'Y1_imp', 'Y2_imp', 'Y3_imp', 'Y4_imp', 'Y5_imp')
n.iters = 3000
n.warmups = 2000
fitStan <- sampling(md, data = stan_dat, iter = n.iters, warmup = n.warmups, 
                    chains = 2, thin=1, init = inits, pars = pars, seed = 2021,
                    control = list(adapt_delta = 0.8, max_treedepth=10))
fname = 'stan_fit/M1_prior2.RData'
save(list = 'fitStan', file = fname)
# load(fname)

A_sample <- extract(fitStan, pars = c('A1', 'A2', 'A3', 'A4', 'A5'))
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

## Filter out longitudinal data with more observations available
nobs = as.vector(table(long$ID))
dat2 = data.frame(ID = unique(long$ID), nobs = nobs)
dat3 = dat2[which(dat2$nobs>=5), ] 
tmp.ID = dat3$ID

## 
Pm = min(round(V/30), 35)
mu_m = colMeans(M.train)
mu_m_smooth = bs.smooth(mu_m, obsgrid, obsgrid, nbasis = Pm)
mu_m_expand = matrix(rep(mu_m_smooth$est.value, N), nrow = N, byrow = T)
M.demean = M.train - mu_m_expand

load('RData/FPCA_Hippo_result.RData')
Lm = length(FPCA.res$values) - L0

phi_m.sign = rep(1, L0)
phi_m.index = rep(1, L0)
psi_m.sign = rep(1, Lm)
psi_m.index = rep(1, Lm)

## Determine the sign of beta_m
t_value = c(0, 0.5, 1, 1.5, 2)
eps = 0.05
cov.mat.full = matrix(NA, nrow = length(t_value)*J, ncol = V)
covariates = matrix(NA, nrow = length(t_value)*J, ncol = L0)
index = 0

for (t0 in t_value){
  long.new = long.i[which(abs(long.i$Time-t0)<eps), ]
  for (j in 1:J){
    index = index + 1
    Y.vec <- as.vector(long.new[, j+2])
    cov.mat.full[index, ] = as.vector(cov(Y.vec, M.train[long.new$ID, ]))
    for (l in 1:L0){
      covariates[index, l] = beta[j]*d0[l]*phi_est[which(tnew==t0), l]
    }
  }
}

H_mat = solve(t(covariates) %*% covariates, t(covariates))
Theta = H_mat %*% cov.mat.full
theta1 = Theta[1, ]
theta2 = Theta[2, ]

phi_m <- matrix(c(theta1, theta2), nrow = V, ncol = L0)
phi_sign <- rep(1, L0)
phi_index <- rep(1, L0)
for (l in 1:L0){
  phi_m[, l] <- sign_eigen(phi_m[, l], phi_index[l], phi_sign[l])
  tmp = sum(phi_m[, l]^2)
  phi_m[, l] <- phi_m[, l]*sqrt(V/tmp)
}

covariates2 = matrix(0, nrow = length(t_value)*J, ncol = V)
index = 0

for (t0 in t_value){
  long.new = long.i[which(abs(long.i$Time-t0)<eps), ]
  for (j in 1:J){
    index = index + 1
    for (l in 1:L0){
      covariates2[index, ] = covariates2[index, ] + 
        beta[j]*d0[l]*phi_est[which(tnew==t0), l]*phi_m[, l]
    } 
  }
}

tmp = as.vector(covariates2)
C_Ym = as.vector(cov.mat.full)
beta_m = unname(coef(lm(C_Ym ~ tmp - 1)))

sign_beta_m = sign(beta_m)

## Estimate eigenfunction of m_i(v) ##
ID1 = tmp.ID
ID2 = tmp.ID
xi1 = summary.xi[ID1, 1:L0]; xi2 = summary.xi[ID2, 1:L0]
M1 = M.demean[ID1, ]; M2 = M.demean[ID2, ]

res = search(M1, M2, xi1, xi2, L0, Lm, V, 
             phi_m.sign, phi_m.index, 
             psi_m.sign, psi_m.index,
             max.iters = 200, eps = 5e-4, 
             beta_m = NULL, Phi_m_est = NULL,
             d0 = d0, sign_beta_m)
beta_m <- res$beta_m
dm <- res$dm
Phi_m_est <- res$Phi_m_est
Psi_m_est <- res$Psi_m_est

## Create m_{ij} matrix ##
m_mat = matrix(NA, nrow = N, ncol = Lm)
f_l = matrix(NA, nrow = L0, ncol = Lm)

for (j in 1:Lm){
  for (i in 1:N){
    m_mat[i, j] <- c_int(M.demean[i, ], Psi_m_est[, j], width = 1/(V-1))
  }
  for (l in 1:L0){
    f_l[l, j] <- c_int(Phi_m_est[, l], Psi_m_est[, j], width = 1/(V-1))
  }
}

## 
md = stan_model('source_code/M1.stan')
stan_dat <- list(n = N, J = J, 
                 nobs = nrow(long), nmiss = Y_missing$len.missing, 
                 id_long = long$ID, 
                 L0 = L0, L1 = L1, Lm = Lm, 
                 P = P, P_surv = ncol(x), 
                 Y1 = Y.list.new[[1]], Y2 = Y.list.new[[2]], Y3 = Y.list.new[[3]],
                 Y4 = Y.list.new[[4]], Y5 = Y.list.new[[5]],
                 miss_index1 = miss.index[[1]], miss_index2 = miss.index[[2]], 
                 miss_index3 = miss.index[[3]], miss_index4 = miss.index[[4]], 
                 miss_index5 = miss.index[[5]],
                 time = long$Time,
                 x = x, surv_time = surv$surv_time, status = surv$status,
                 b = bs.mat, phi = phi, psi = psi, 
                 m_mat = m_mat, f_l = f_l, beta_m = beta_m, 
                 Ltau = n.tau, tau = tau[1:n.tau], 
                 h_grid = h.grid, h_index = h.index)
rnd <- 1/3
inits1 <- list(A1 = g.inits(mu_est[[1]]$coefficient), A2 = g.inits(mu_est[[2]]$coefficient),
               A3 = g.inits(mu_est[[3]]$coefficient), A4 = g.inits(mu_est[[4]]$coefficient), 
               A5 = g.inits(mu_est[[5]]$coefficient), 
               sqrt_d0 = as.array(g.inits.sd(sqrt(d0))), 
               sqrt_d1 = as.array(g.inits.sd(sqrt(d1))),
               xi = t(summary.xi[, 1:L0]), 
               beta = g.inits(beta[2:J]), 
               sigma = g.inits.sd(sigma),
               logh0 = g.inits(logh0_hat), 
               gamma_x = as.array(g.inits(gamma_x_hat)), 
               gamma0 = as.array(rnorm(L0, 0, 1/3)),
               gamma11 = as.array(rnorm(L1, 0, 1/3)), gamma12 = as.array(rnorm(L1, 0, 1/3)),
               gamma13 = as.array(rnorm(L1, 0, 1/3)), gamma14 = as.array(rnorm(L1, 0, 1/3)),
               gamma15 = as.array(rnorm(L1, 0, 1/3)), 
               gamma_m = as.array(rnorm(Lm, 0, 1/3)), 
               Y1_imp = mean_Y[[1]]+rnd, Y2_imp = mean_Y[[2]]+rnd, Y3_imp = mean_Y[[3]]+rnd, 
               Y4_imp = mean_Y[[4]]+rnd, Y5_imp = mean_Y[[5]]+rnd)
inits2 <- list(A1 = g.inits(mu_est[[1]]$coefficient), A2 = g.inits(mu_est[[2]]$coefficient),
               A3 = g.inits(mu_est[[3]]$coefficient), A4 = g.inits(mu_est[[4]]$coefficient), 
               A5 = g.inits(mu_est[[5]]$coefficient), 
               sqrt_d0 = as.array(g.inits.sd(sqrt(d0))), 
               sqrt_d1 = as.array(g.inits.sd(sqrt(d1))),
               xi = t(summary.xi[, 1:L0]), 
               beta = g.inits(beta[2:J]), 
               sigma = g.inits.sd(sigma),
               logh0 = g.inits(logh0_hat), 
               gamma_x = as.array(g.inits(gamma_x_hat)), 
               gamma0 = as.array(rnorm(L0, 0, 1/3)),
               gamma11 = as.array(rnorm(L1, 0, 1/3)), gamma12 = as.array(rnorm(L1, 0, 1/3)),
               gamma13 = as.array(rnorm(L1, 0, 1/3)), gamma14 = as.array(rnorm(L1, 0, 1/3)),
               gamma15 = as.array(rnorm(L1, 0, 1/3)), 
               gamma_m = as.array(rnorm(Lm, 0, 1/3)), 
               Y1_imp = mean_Y[[1]]+rnd, Y2_imp = mean_Y[[2]]+rnd, Y3_imp = mean_Y[[3]]+rnd, 
               Y4_imp = mean_Y[[4]]+rnd, Y5_imp = mean_Y[[5]]+rnd)
inits <- list(c1 = inits1, c2 = inits2)
pars <- c('A1', 'A2', 'A3', 'A4', 'A5', 
          'd0', 'd1', 'dm', 
          'beta', 'sigma',
          'logh0', 'gamma_x', 'gamma0', 
          'gamma11', 'gamma12', 'gamma13',  'gamma14', 'gamma15', 'gamma_m', 
          'xi', 'zeta_1', 'zeta_2', 'zeta_3', 'zeta_4', 'zeta_5', 'LL', 
          'Y1_imp', 'Y2_imp', 'Y3_imp', 'Y4_imp', 'Y5_imp')
n.iters = 5000
n.warmups = 4000
fitStan <- sampling(md, data = stan_dat, iter = n.iters, warmup = n.warmups, 
                    chains = 2, thin=1, init = inits, pars = pars, seed = 2021,
                    control = list(adapt_delta = 0.8, max_treedepth=10))
fname <- paste0('stan_fit/', 'M1_posterior.RData')
save(list = 'fitStan', file = fname)

summ = summary(fitStan)$summary[1:(J*P + L0 + L1 + Lm + J-1 + J + n.tau + ncol(x) + L0 + L1*J + Lm), c(1, 3, 4, 8, 9, 10)]
write.csv(summ, file = 'RData/summary_posterior_M1.csv')

fname <- paste0('RData/summary_posterior_M1.RData')
save(list = c('beta_m', 'phi_est', 'psi_est', 'Phi_m_est', 'Psi_m_est'), file = fname)

## Compute log likelihood ##
Q <- 2000 ## 2000 samples
A_sample <- extract(fitStan, pars = c('A1', 'A2', 'A3', 'A4', 'A5'))
d0_sample <- extract(fitStan, pars = 'd0')$d0
d1_sample <- extract(fitStan, pars = 'd1')$d1
dm_sample <- extract(fitStan, pars = 'dm')$dm
beta_sample <- cbind(rep(1, Q), extract(fitStan, pars = 'beta')$beta)
sigma_sample <- extract(fitStan, par = 'sigma')$sigma
logh0_sample <- extract(fitStan, 'logh0')$logh0
gamma_x_sample <- extract(fitStan, pars = 'gamma_x')$gamma_x
gamma0_sample <- extract(fitStan, pars = 'gamma0')$gamma0
gamma11_sample <- extract(fitStan, pars = 'gamma11')$gamma11
gamma12_sample <- extract(fitStan, pars = 'gamma12')$gamma12
gamma13_sample <- extract(fitStan, pars = 'gamma13')$gamma13
gamma14_sample <- extract(fitStan, pars = 'gamma14')$gamma14
gamma15_sample <- extract(fitStan, pars = 'gamma15')$gamma15
gamma_m_sample <- extract(fitStan, pars = 'gamma_m')$gamma_m
xi_sample <- extract(fitStan, pars = 'xi')$xi ## Q*L0*N
zeta1_sample <- extract(fitStan, pars = 'zeta_1')$zeta_1
zeta2_sample <- extract(fitStan, pars = 'zeta_2')$zeta_2
zeta3_sample <- extract(fitStan, pars = 'zeta_3')$zeta_3
zeta4_sample <- extract(fitStan, pars = 'zeta_4')$zeta_4
zeta5_sample <- extract(fitStan, pars = 'zeta_5')$zeta_5
LL_sample <- extract(fitStan, pars = 'LL')$LL
Y_imp_sample <- extract(fitStan, pars = c('Y1_imp', 'Y2_imp', 'Y3_imp', 'Y4_imp', 'Y5_imp'))

ll_surv <- LL_sample

## Calculation of log likelihood ##
nobs <- nrow(long); ID <- long$ID
ll_long <- array(NA, dim = c(Q, nobs, J)) ## Q*nobs*J
ll_long_full <- array(0, dim = c(Q, N)) ## Q*N
for (q in 1:Q){
  for (i in 1:nobs){
    Y.list.tmp <- Y.list.o
    for (j in 1:J){
      if (i %in% miss.index[[j]]){
        tmp.index <- which(i==miss.index[[j]])
        Y.list.tmp[[j]][i] <- Y_imp_sample[[j]][q, tmp.index]
      }
    }
    
    u <- sum(xi_sample[q, , ID[i]]*phi[i, ])
    w1 <- sum(zeta1_sample[q, , ID[i]]*psi[i, ])
    w2 <- sum(zeta2_sample[q, , ID[i]]*psi[i, ])
    w3 <- sum(zeta3_sample[q, , ID[i]]*psi[i, ])
    w4 <- sum(zeta4_sample[q, , ID[i]]*psi[i, ])
    w5 <- sum(zeta5_sample[q, , ID[i]]*psi[i, ])
    w <- list(w1 = w1, w2 = w2, w3 = w3, w4 = w4, w5 = w5)
    for (j in 1:J){
      mu_j <- sum(A_sample[[j]][q, ]*bs.mat[i, ]) + beta_sample[q, j]*(u + w[[j]])
      ll_long[q, i, j] <- dsn(Y.list.tmp[[j]][i], mu_j, 
                              sigma_sample[q, j], 0,  log = T)
    }
    
    ll_long_full[q, ID[i]] <-  ll_long_full[q, ID[i]] + sum(ll_long[q, i, ])
  }
}

Dbar.long <- sum(-2*ll_long_full)/Q
Dbar.surv <- sum(-2*ll_surv)/Q
Dbar <- Dbar.long + Dbar.surv

np <- nrow(summ) - Lm    ## dm is not a parameter
EAIC <- Dbar + 2*np
EBIC <- Dbar + log(N)*np

## Compute LOOIC, WAIC ##
ll_full <- ll_long_full + ll_surv
rel_n_eff <- relative_eff(exp(ll_full), chain_id = rep(1:2, each = 1000))
looic <- loo(ll_full, r_eff = rel_n_eff, cores = 4)$estimates[3, 1]
waic <- waic(ll_full)$estimates[3, 1]

l <- list(Dbar.long = Dbar.long, Dbar.surv = Dbar.surv, Dbar = Dbar, 
          np = np, EAIC = EAIC, EBIC = EBIC,
          looic = looic, waic = waic)
save(list = 'l', file = 'RData/summary_posterior_M1_MP.RData')
