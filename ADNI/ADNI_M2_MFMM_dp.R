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
library(statmod)
library(tdROC)
library(ipred)
rstan_options(auto_write = TRUE)

# seed <- 1
set.seed(2011*seed)

source('source_code/functions.R')
source('source_code/sFPCA.R')

long <- read.csv('dataset/ADNI_long_2.csv')
surv <- read.csv('dataset/ADNI_surv.csv')

N <- nrow(surv)
J <- 5 ## number of longitudinal outcomes
P <- 9 ## number of basis functions for B-spline
L0 <- 2 ## npc for U_i(t)
L1 <- 1 ## npc for W_ij(t)

## Training and testing set
N.train <- round(N*0.75); N.test <- N - N.train; ID <- 1:N
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

## Refine the time grid ##
tnew <- (0:101)/10
tg <- length(tnew)

fit.cox <- coxph(Surv(surv_time, status) ~ Age + Gender + Education + APOE4,
                 data = surv.train, method = 'breslow')
gamma_x_hat <- unname(fit.cox$coefficients)
baseline.cumulative.hazard <- basehaz(fit.cox)
x <- as.matrix(cbind(surv.train$Age, surv.train$Gender, surv.train$Education, surv.train$APOE4))

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

## Use mfaces package to calculate covariance matrices ##
dat.mface <- list('y1' = data.frame('subj' = long.train$ID, 'argvals' = long.train$Time, 'y' = Y.list.i[[1]]),
                  'y2' = data.frame('subj' = long.train$ID, 'argvals' = long.train$Time, 'y' = Y.list.i[[2]]),
                  'y3' = data.frame('subj' = long.train$ID, 'argvals' = long.train$Time, 'y' = Y.list.i[[3]]),
                  'y4' = data.frame('subj' = long.train$ID, 'argvals' = long.train$Time, 'y' = Y.list.i[[4]]),
                  'y5' = data.frame('subj' = long.train$ID, 'argvals' = long.train$Time, 'y' = Y.list.i[[5]]))
fit.mface <- mface.sparse(dat.mface, argvals.new = tnew, knots = 6, newdata = dat.mface)
# save(list = 'fit.mface', file = 'RData/fit_mface_dp1.RData')

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
phi <- matrix(NA, nrow(long.train), ncol = L0)
phi_surv <- matrix(NA, N.train, ncol = L0) ## phi_l(T_i) used in evaluation of h_i(T_i) 
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

## Create Gaussian Quadrature points and weights at interval (a_{ie}, b_{ie})
create_gq <- function(weights, nodes, time0, time1){
  const <- rep(NA, length(time0)) ## const_{i} = (b_{i} - a_{i})/2
  a <- b <- rep(NA, length(time0))
  
  phi1_interval <- phi2_interval <- array(NA, dim = c(length(time0), G))
  psi1_interval <- array(NA, dim = c(length(time0), G))
  for (i in 1:length(time0)){
    a[i] <- time0[i]
    b[i] <- time1[i]
    const[i] <- (b[i] - a[i])/2
    x_star <- (b[i] - a[i])/2*nodes + (b[i] + a[i])/2
    
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
weights <- gauss.quad(G, kind = 'legendre')$weights ## weights
nodes <- gauss.quad(G, kind = 'legendre')$nodes  ## nodes
const <- matrix(NA, nrow = N.train, ncol = n.tau) ## const_{ie} = (b_{ie} - a_{ie})/2

## We evaluate phi_l(x_star) and psi_l(x_star) 
## where x_star = (b[i, e] - a[i, e])/2*x_g + (b[i, e] + a[i, e])/2
phi1_interval <- phi2_interval <- array(NA, dim = c(N.train, n.tau, G))
psi1_interval <- array(NA, dim = c(N.train, n.tau, G))

for (e in 1:n.tau){
  l <- create_gq(weights, nodes, h.grid[, e], h.grid[, e+1])
  const[, e] <- l$const
  phi1_interval[, e, ] <- l$phi1_interval
  phi2_interval[, e, ] <- l$phi2_interval
  psi1_interval[, e, ] <- l$psi1_interval
}

## Stan fit ##
md = stan_model('source_code/M2_prior1.stan')
stan_dat <- list(n = N.train, J = J, 
                 nobs = nrow(long.train), nmiss = Y_missing$len.missing, 
                 id_long = long.train$ID, 
                 L0 = L0, L1 = L1, P = P, P_surv = ncol(x), 
                 Y1 = Y.list.new[[1]], Y2 = Y.list.new[[2]], Y3 = Y.list.new[[3]],
                 Y4 = Y.list.new[[4]], Y5 = Y.list.new[[5]],
                 miss_index1 = miss.index[[1]], miss_index2 = miss.index[[2]], 
                 miss_index3 = miss.index[[3]], miss_index4 = miss.index[[4]], 
                 miss_index5 = miss.index[[5]],
                 time = long.train$Time, x = x, 
                 surv_time = surv.train$surv_time, status = surv.train$status,
                 b = bs.mat, phi = phi, psi = psi,
                 phi_surv = phi_surv, psi_surv = psi_surv,
                 Ltau = n.tau, tau = tau[1:n.tau], 
                 h_grid = h.grid, h_index = h.index,
                 G = G, w = weights, constant = const, 
                 phi1_interval = phi1_interval, phi2_interval = phi2_interval,
                 psi1_interval = psi1_interval)

set.seed(2021*seed)
rnd <- 1/3
inits1 <- list(A1 = g.inits(mu_est[[1]]$coefficient), A2 = g.inits(mu_est[[2]]$coefficient),
               A3 = g.inits(mu_est[[3]]$coefficient), A4 = g.inits(mu_est[[4]]$coefficient), 
               A5 = g.inits(mu_est[[5]]$coefficient), 
               sqrt_d0 = as.array(g.inits.sd(sqrt(d0_est))), 
               sqrt_d1 = as.array(g.inits.sd(sqrt(d1_est))),
               beta = g.inits(beta_hat[2:J]), 
               sigma = g.inits.sd(sigma_hat),
               logh0 = g.inits(logh0_hat), 
               gamma_x = as.array(g.inits(gamma_x_hat)), 
               gamma0 = rnorm(1, 0, 1/3),
               gamma1 = as.array(rnorm(J, 0, 1/3)),  
               Y1_imp = mean_Y[[1]]+rnd, Y2_imp = mean_Y[[2]]+rnd, Y3_imp = mean_Y[[3]]+rnd, 
               Y4_imp = mean_Y[[4]]+rnd, Y5_imp = mean_Y[[5]]+rnd)
inits2 <- list(A1 = g.inits(mu_est[[1]]$coefficient), A2 = g.inits(mu_est[[2]]$coefficient),
               A3 = g.inits(mu_est[[3]]$coefficient), A4 = g.inits(mu_est[[4]]$coefficient), 
               A5 = g.inits(mu_est[[5]]$coefficient), 
               sqrt_d0 = as.array(g.inits.sd(sqrt(d0_est))), 
               sqrt_d1 = as.array(g.inits.sd(sqrt(d1_est))),
               beta = g.inits(beta_hat[2:J]), 
               sigma = g.inits.sd(sigma_hat),
               logh0 = g.inits(logh0_hat), 
               gamma_x = as.array(g.inits(gamma_x_hat)), 
               gamma0 = rnorm(1, 0, 1/3),
               gamma1 = as.array(rnorm(J, 0, 1/3)),  
               Y1_imp = mean_Y[[1]]-rnd, Y2_imp = mean_Y[[2]]-rnd, Y3_imp = mean_Y[[3]]-rnd, 
               Y4_imp = mean_Y[[4]]-rnd, Y5_imp = mean_Y[[5]]-rnd)

inits <- list(c1 = inits1, c2 = inits2)
pars <- c('A1', 'A2', 'A3', 'A4', 'A5', 
          'd0', 'd1', 'beta', 
          'sigma', 
          'logh0', 'gamma_x', 'gamma0', 
          'gamma1',  'xi')
n.iters = 3000
n.warmups = 2000
fitStan <- sampling(md, data = stan_dat, iter = n.iters, warmup = n.warmups, 
                    chains = 2, thin=1, init = inits, pars = pars, seed = 2021,
                    control = list(adapt_delta = 0.8, max_treedepth=10))
fname = paste0('stan_fit/M2_MFMM_dp/', seed, '_prior1.RData')
save(list = 'fitStan', file = fname)
# load(fname)

## Extract samples ##
Q <- 2000 ## 2000 samples
A_sample <- extract(fitStan, pars = c('A1', 'A2', 'A3', 'A4', 'A5'))
d0_sample <- extract(fitStan, pars = 'd0')$d0
d1_sample <- extract(fitStan, pars = 'd1')$d1
beta_sample <- extract(fitStan, pars = 'beta')$beta
sigma_sample <- extract(fitStan, pars = 'sigma')$sigma
logh0_sample <- extract(fitStan, pars = 'logh0')$logh0
gamma_x_sample <- extract(fitStan, pars = 'gamma_x')$gamma_x
gamma0_sample <- extract(fitStan, pars = 'gamma0')$gamma0
gamma1_sample <- extract(fitStan, pars = 'gamma1')$gamma1

## Compute posterior mean
A <- lapply(1:J, FUN = function(x) colMeans(A_sample[[x]]))
d0 <- apply(d0_sample, 2, 'mean')
d1 <- apply(d1_sample, 2, 'mean')
beta <- apply(beta_sample, 2, 'mean')
sigma <- sapply(1:J, FUN = function(x) mean(sigma_sample[, x]))
logh0 <- apply(logh0_sample, 2, 'mean')
gamma_x <- apply(gamma_x_sample, 2, 'mean')
gamma0 <- mean(gamma0_sample)
gamma1 <- apply(gamma1_sample, 2, 'mean')

## Dynamic prediction settings
md_dp = stan_model('source_code/M2_MFMM_dp.stan')
pars_dp = c('xi', 'zeta_1', 'zeta_2', 'zeta_3', 'zeta_4', 'zeta_5')
starting.time <- c(2, 2.5, 3, 3.5, 4)
delta.time <- seq(0.5, 1.5, by = 0.01)
min.log <- log(.Machine$double.xmin) + 10
AUC <- BS <- matrix(NA, length(starting.time), length(delta.time))

## compute estimated survival probability at time t: exp(-cumulative hazard)
H_est <- function(phi1_interval, phi2_interval, psi1_interval, constant, 
                  x, weights, h_grid, xi, zeta1, zeta2, zeta3, zeta4, zeta5, t){
  
  H <- rep(NA, n.tau)
  for (k in 1:n.tau){
    U = xi[1]*phi1_interval[k, ] + xi[2]*phi2_interval[k, ]
    W1 = zeta1*psi1_interval[k, ]
    W2 = zeta2*psi1_interval[k, ]
    W3 = zeta3*psi1_interval[k, ]
    W4 = zeta4*psi1_interval[k, ]
    W5 = zeta5*psi1_interval[k, ]
    f = exp(gamma0*U + gamma1[1]*W1 + gamma1[2]*W2 + gamma1[3]*W3 + gamma1[4]*W4 + gamma1[5]*W5);
    H[k] = exp(sum(x*gamma_x) + logh0[k])*constant[k]*sum(f*weights);
  }
  return(sum(H))
}

for (Tstart in starting.time){
  
  ## filter out subjects with observed survival time greater than starting time
  surv.test2 <- surv.test[surv.test$surv_time>Tstart, ]
  long.test2 <- long.test[long.test$ID %in% surv.test2$ID, ]
  surv.test2$newID <- 1:nrow(surv.test2)
  long.test2$newID <- rep(1:nrow(surv.test2), table(long.test2$ID))
  x.test <- as.matrix(cbind(surv.test2$Age, surv.test2$Gender, surv.test2$Education, surv.test2$APOE4))
  
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
  Y.list.test <- list(y1 = long.test.prior$y1, y2 = long.test.prior$y2, y3 = long.test.prior$y3, 
                      y4 = long.test.prior$y4, y5 = long.test.prior$y5)
  Y_missing.test <- missing_determine(value = Y.list.test)
  
  h_grid_T0 <- rep(NA, n.tau)
  for (k in 1:n.tau){
    h_grid_T0[k] <- ifelse(Tstart>=tau[k], min(Tstart, tau[k+1]), Tstart)
  }
  h_grid_T0 <- c(0, h_grid_T0)
  
  const_T0 <- rep(NA, n.tau) 
  phi1_interval_T0 <- phi2_interval_T0 <- array(NA, dim = c(n.tau, G))
  psi1_interval_T0 <- array(NA, dim = c(n.tau, G))
  
  for (e in 1:n.tau){
    l <- create_gq(weights, nodes, h_grid_T0[e], h_grid_T0[e+1])
    const_T0[e] <- l$const
    phi1_interval_T0[e, ] <- l$phi1_interval
    phi2_interval_T0[e, ] <- l$phi2_interval
    psi1_interval_T0[e, ] <- l$psi1_interval
  }
  
  # NUTS
  stan_dat_dp <- list(n = length(tmp.ID.test), J = J, 
                      nobs = nrow(long.test.prior),
                      id_long = long.test.prior$newID,
                      L0 = L0, L1 = L1, 
                      P = P, P_surv = ncol(x.test), 
                      Y1 = Y_missing.test$new.value[[1]], 
                      Y2 = Y_missing.test$new.value[[2]], 
                      Y3 = Y_missing.test$new.value[[3]],
                      Y4 = Y_missing.test$new.value[[4]],
                      Y5 = Y_missing.test$new.value[[5]],
                      time = long.test.prior$Time, x = x.test,
                      surv_time = rep(Tstart, length(tmp.ID.test)), 
                      b = bs.mat.test, phi = phi_test, psi = psi_test,
                      Ltau = n.tau, tau = tau[1:n.tau], 
                      h_grid = h_grid_T0, 
                      G = G, w = weights, constant = const_T0, 
                      phi1_interval = phi1_interval_T0, phi2_interval = phi2_interval_T0,
                      psi1_interval = psi1_interval_T0,
                      A1 = A[[1]], A2 = A[[2]], A3 = A[[3]], A4 = A[[4]], A5 = A[[5]],
                      sqrt_d0 = as.array(sqrt(d0)), sqrt_d1 = as.array(sqrt(d1)), 
                      beta = beta, sigma = sigma,
                      logh0 = logh0, gamma_x = as.array(gamma_x), 
                      gamma0 = gamma0, 
                      gamma1 = gamma1)
  inits1_dp <- inits2_dp <- list(xi = matrix(0, L0, length(tmp.ID.test)), 
                                 zeta1 = matrix(0, L1, length(tmp.ID.test)), 
                                 zeta2 = matrix(0, L1, length(tmp.ID.test)), 
                                 zeta3 = matrix(0, L1, length(tmp.ID.test)), 
                                 zeta4 = matrix(0, L1, length(tmp.ID.test)), 
                                 zeta5 = matrix(0, L1, length(tmp.ID.test)))
  inits_dp <- list(c1 = inits1_dp, c2 = inits2_dp)
  fitStan_dp <- sampling(md_dp, data = stan_dat_dp, iter = 2000, warmup = 1000, 
                         chains = 2, thin=1, init = inits_dp, pars = pars_dp, seed = 123,
                         control = list(adapt_delta = 0.8, max_treedepth=10))
  xi_sample <- extract(fitStan_dp, pars = 'xi')$xi ## extract xi_{il}
  zeta1_sample <- extract(fitStan_dp, pars = 'zeta_1')$zeta_1 ## extract zeta1_{il}
  zeta2_sample <- extract(fitStan_dp, pars = 'zeta_2')$zeta_2 ## extract zeta2_{il}
  zeta3_sample <- extract(fitStan_dp, pars = 'zeta_3')$zeta_3
  zeta4_sample <- extract(fitStan_dp, pars = 'zeta_4')$zeta_4
  zeta5_sample <- extract(fitStan_dp, pars = 'zeta_5')$zeta_5
  
  ## calculate conditional survival probabilities at T+\delta_T
  for (Tdelta in delta.time){
    
    h_grid_T1 <- rep(NA, n.tau)
    for (k in 1:n.tau){
      h_grid_T1[k] <- ifelse(Tstart+Tdelta>=tau[k], min(Tstart+Tdelta, tau[k+1]), Tstart+Tdelta)
    }
    h_grid_T1 <- c(0, h_grid_T1)
    
    const_T1 <- rep(NA, n.tau) 
    phi1_interval_T1 <- phi2_interval_T1 <- array(NA, dim = c(n.tau, G))
    psi1_interval_T1 <- array(NA, dim = c(n.tau, G))
    
    for (e in 1:n.tau){
      l <- create_gq(weights, nodes, h_grid_T1[e], h_grid_T1[e+1])
      const_T1[e] <- l$const
      phi1_interval_T1[e, ] <- l$phi1_interval
      phi2_interval_T1[e, ] <- l$phi2_interval
      psi1_interval_T1[e, ] <- l$psi1_interval
    }
    
    tmp.surv.predict <- rep(NA, length(tmp.ID.test))
    
    for (i in 1:length(tmp.ID.test)){
      xi_i <- xi_sample[, , i] 
      zeta1_i <- zeta1_sample[, , i]
      zeta2_i <- zeta2_sample[, , i]
      zeta3_i <- zeta3_sample[, , i]
      zeta4_i <- zeta4_sample[, , i]
      zeta5_i <- zeta5_sample[, , i]
      Hi_est_t0 <- Hi_est_t1 <- rep(NA, Q)
      for (q in 1:Q){
        Hi_est_t0[q] <- H_est(phi1_interval_T0, phi2_interval_T0, psi1_interval_T0, 
                              const_T0, x.test[i, ], weights, h_grid_T0, 
                              xi_i[q, ], zeta1_i[q], zeta2_i[q], zeta3_i[q], 
                              zeta4_i[q], zeta5_i[q], 
                              t = Tstart)
        Hi_est_t1[q] <- H_est(phi1_interval_T1, phi2_interval_T1, psi1_interval_T1, 
                              const_T1, x.test[i, ], weights, h_grid_T1, 
                              xi_i[q, ], zeta1_i[q], zeta2_i[q], zeta3_i[q], 
                              zeta4_i[q], zeta5_i[q],
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

fname <- paste0('result/M2_MFMM_dp/', seed, '.RData')
save(list = c('AUC', 'BS'), file = fname)
