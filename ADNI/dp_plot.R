## Dynamic Prediction plot ##
options(mc.cores = parallel::detectCores())
library(tidyverse)
library(data.table)
library(splines)
library(survival)
library(rstan)
library(Matrix)
library(refund)
library(orthogonalsplinebasis)
library(mfaces)
library(statmod)
library(loo)
library(Amelia)
library(sn)
library(dlm)
library(tdROC)
library(ipred)
rstan_options(auto_write = TRUE)

# seed <- 1
set.seed(2031*seed)

source('source_code/functions.R')
source('source_code/sFPCA.R')
source('source_code/search.R')

long <- read.csv('dataset/ADNI_long_2.csv')
surv <- read.csv('dataset/ADNI_surv.csv')
MRI.dat <- fread("dataset/Surv.mri.dat.csv") %>% as.matrix()

N <- nrow(surv)
J <- 5 ## number of longitudinal outcomes
V = ncol(MRI.dat)
obsgrid = (1:V)/V
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

M.train = MRI.dat[ID.train, ]
M.test = MRI.dat[ID.test, ]

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
md = stan_model('source_code/M2_prior2.stan')
stan_dat <- list(n = N.train, J = J, 
                 nobs = nrow(long.train), nmiss = Y_missing$len.missing, 
                 id_long = long.train$ID, 
                 L0 = L0, L1 = L1, P = P, 
                 Y1 = Y.list.new[[1]], Y2 = Y.list.new[[2]], Y3 = Y.list.new[[3]],
                 Y4 = Y.list.new[[4]], Y5 = Y.list.new[[5]],
                 miss_index1 = miss.index[[1]], miss_index2 = miss.index[[2]], 
                 miss_index3 = miss.index[[3]], miss_index4 = miss.index[[4]], 
                 miss_index5 = miss.index[[5]],
                 time = long.train$Time, 
                 b = bs.mat, phi = phi, psi = psi)

set.seed(2021*seed)
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
          'sigma', 'xi')
n.iters = 3000
n.warmups = 2000
fitStan <- sampling(md, data = stan_dat, iter = n.iters, warmup = n.warmups, 
                    chains = 2, thin=1, init = inits, pars = pars, seed = 2021,
                    control = list(adapt_delta = 0.8, max_treedepth=10))
fname = paste0('stan_fit/M2_whole_dp/', seed, '_prior2.RData')
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
nobs = as.vector(table(long.train$ID))
dat2 = data.frame(ID = unique(long.train$ID), nobs = nobs)
dat3 = dat2[which(dat2$nobs>=5), ] 
tmp.ID = dat3$ID

## 
Pm = min(round(V/30), 35)
mu_m = colMeans(M.train)
mu_m_smooth = bs.smooth(mu_m, obsgrid, obsgrid, nbasis = Pm)
mu_m_expand = matrix(rep(mu_m_smooth$est.value, N.train), nrow = N.train, byrow = T)
M.demean = M.train - mu_m_expand

Lm = 5

phi_m.sign = rep(1, L0)
phi_m.index = rep(1, L0)
psi_m.sign = rep(1, Lm)
psi_m.index = rep(1, Lm)

sign_beta_m = -1

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
                 P = P, P_surv = ncol(x), 
                 Y1 = Y.list.new[[1]], Y2 = Y.list.new[[2]], Y3 = Y.list.new[[3]],
                 Y4 = Y.list.new[[4]], Y5 = Y.list.new[[5]],
                 miss_index1 = miss.index[[1]], miss_index2 = miss.index[[2]], 
                 miss_index3 = miss.index[[3]], miss_index4 = miss.index[[4]], 
                 miss_index5 = miss.index[[5]],
                 time = long.train$Time,
                 x = x, surv_time = surv.train$surv_time, status = surv.train$status,
                 b = bs.mat, phi = phi, psi = psi, 
                 phi_surv = phi_surv, psi_surv = psi_surv,
                 m_mat = m_mat, f_l = f_l, beta_m = beta_m, 
                 Ltau = n.tau, tau = tau[1:n.tau], 
                 h_grid = h.grid, h_index = h.index,
                 G = G, w = weights, constant = const, 
                 phi1_interval = phi1_interval, phi2_interval = phi2_interval,
                 psi1_interval = psi1_interval)
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
               gamma0 = rnorm(1, 0, 1/3),
               gamma1 = as.array(rnorm(J, 0, 1/3)),  
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
               gamma0 = rnorm(1, 0, 1/3),
               gamma1 = as.array(rnorm(J, 0, 1/3)),  
               gamma_m = as.array(rnorm(Lm, 0, 1/3)),  
               Y1_imp = mean_Y[[1]]+rnd, Y2_imp = mean_Y[[2]]+rnd, Y3_imp = mean_Y[[3]]+rnd, 
               Y4_imp = mean_Y[[4]]+rnd, Y5_imp = mean_Y[[5]]+rnd)
inits <- list(c1 = inits1, c2 = inits2)
pars <- c('A1', 'A2', 'A3', 'A4', 'A5', 
          'd0', 'd1', 'dm', 
          'beta', 'sigma',
          'logh0', 'gamma_x', 'gamma0', 
          'gamma1', 'gamma_m')
n.iters = 5000
n.warmups = 4000
fitStan <- sampling(md, data = stan_dat, iter = n.iters, warmup = n.warmups, 
                    chains = 2, thin=1, init = inits, pars = pars, seed = 2021,
                    control = list(adapt_delta = 0.8, max_treedepth=10))
fname = paste0('stan_fit/M2_whole_dp/', seed, '_posterior.RData')
save(list = 'fitStan', file = fname)

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
gamma_m_sample <- extract(fitStan, pars = 'gamma_m')$gamma_m

## Compute posterior mean
A <- lapply(1:J, FUN = function(x) colMeans(A_sample[[x]]))
d0 <- apply(d0_sample, 2, 'mean')
d1 <- apply(d1_sample, 2, 'mean')
beta <- apply(beta_sample, 2, 'mean')
beta <- c(1, beta)
sigma <- sapply(1:J, FUN = function(x) mean(sigma_sample[, x]))
logh0 <- apply(logh0_sample, 2, 'mean')
gamma_x <- apply(gamma_x_sample, 2, 'mean')
gamma0 <- mean(gamma0_sample)
gamma1 <- apply(gamma1_sample, 2, 'mean')
gamma_m <- apply(gamma_m_sample, 2, 'mean')

md_dp = stan_model('source_code/M2_dp.stan')
pars_dp = c('xi', 'zeta_1', 'zeta_2', 'zeta_3', 'zeta_4', 'zeta_5', 'xi_m')
starting.time <- c(2, 3)

## compute estimated survival probability at time t: exp(-cumulative hazard)
H_est <- function(phi1_interval, phi2_interval, psi1_interval, constant, 
                  x, weights, h_grid, xi, zeta1, zeta2, zeta3, zeta4, zeta5, xi_m, t){
  
  H <- rep(NA, n.tau)
  for (k in 1:n.tau){
    U = xi[1]*phi1_interval[k, ] + xi[2]*phi2_interval[k, ]
    W1 = zeta1*psi1_interval[k, ]
    W2 = zeta2*psi1_interval[k, ]
    W3 = zeta3*psi1_interval[k, ]
    W4 = zeta4*psi1_interval[k, ]
    W5 = zeta5*psi1_interval[k, ]
    f = exp(gamma0*U + gamma1[1]*W1 + gamma1[2]*W2 + gamma1[3]*W3 + gamma1[4]*W4 + gamma1[5]*W5);
    H[k] = exp(sum(x*gamma_x) + logh0[k] + sum(xi_m*gamma_m))*constant[k]*sum(f*weights);
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
  
  M.test2 <- M.test[surv.test2$ID - N.train, ]
  mu_m_expand = matrix(rep(mu_m_smooth$est.value, nrow(surv.test2)), nrow = nrow(surv.test2), byrow = T)
  M.demean = M.test2 - mu_m_expand
  
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
                      L0 = L0, L1 = L1, Lm = Lm, 
                      P = P, P_surv = ncol(x.test), 
                      Y1 = Y_missing.test$new.value[[1]], 
                      Y2 = Y_missing.test$new.value[[2]], 
                      Y3 = Y_missing.test$new.value[[3]],
                      Y4 = Y_missing.test$new.value[[4]],
                      Y5 = Y_missing.test$new.value[[5]],
                      time = long.test.prior$Time, x = x.test,
                      surv_time = rep(Tstart, length(tmp.ID.test)), 
                      b = bs.mat.test, phi = phi_test, psi = psi_test,
                      m_mat = m_mat, f_l = f_l, beta_m = beta_m, 
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
                      gamma1 = gamma1, 
                      gamma_m = as.array(gamma_m))
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
  xi_m_sample <- extract(fitStan_dp, pars = 'xi_m')$xi_m ## extract xi_{il}
  fname <- paste0('stan_fit/M2_whole_dp/dp1_', which(Tstart==starting.time), '.RData')
  save(list = c('xi_sample', 'zeta1_sample', 'zeta2_sample', 
                'zeta3_sample', 'zeta4_sample', 'zeta5_sample', 'xi_m_sample'), file = fname)
  rm(fitStan_dp)
}

# ID.list = c(559, 561, 566, 571, 578, 586, 603, 
#             608, 611, 635, 636, 642, 644, 647, 
#             655, 658, 659, 664, 667, 669, 677,
#             679, 690, 699, 701, 704, 708, 709,
#             712, 719, 727, 742)
ID.list = c(572, 593, 596, 599, 625, 629, 637, 
            645, 649, 657, 695, 711, 729, 731)

for (tmp.ID.new in ID.list){

for (Tstart in starting.time){
  
  delta.time <- seq(0, 8, by = 0.1)
  delta.time <- delta.time[which(delta.time+Tstart<=10.1)]
  
  setEPS()
  fname <- paste0('plot/dp2/', tmp.ID.new, '_', which(Tstart==starting.time), '.eps')
  postscript(fname, height = 5)
  
  par(mfrow = c(2, 3),     
      oma = c(0, 0, 0, 0),
      mar = c(3, 4, 2, 2))
  surv.test2 <- surv.test[surv.test$surv_time>Tstart, ]
  long.test2 <- long.test[long.test$ID %in% surv.test2$ID, ]
  surv.test2$newID <- 1:nrow(surv.test2)
  long.test2$newID <- rep(1:nrow(surv.test2), table(long.test2$ID))
  
  x.test <- as.matrix(cbind(surv.test2$Age, surv.test2$Gender, surv.test2$Education, surv.test2$APOE4))
  
  M.test2 <- M.test[surv.test2$ID - N.train, ]
  mu_m_expand = matrix(rep(mu_m_smooth$est.value, nrow(surv.test2)), nrow = nrow(surv.test2), byrow = T)
  M.demean = M.test2 - mu_m_expand
  
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
  
  t.new <- Tstart + delta.time
  bs.mat.new <- create_bs(time.grid = tnew, pred.time = t.new, nbasis = P)
  
  phi.new <- matrix(NA, length(t.new), ncol = L0)
  for (l in 1:L0){
    phi.smooth <- bs.smooth(phi_est[, l], tnew, argvals.new = t.new, nbasis = P)
    phi.new[, l] <- phi.smooth$est.value
  }
  
  psi.new <- matrix(NA, length(t.new), ncol = L1)
  for (l in 1:L1){
    psi.smooth <- bs.smooth(psi_est[, l], tnew, argvals.new = t.new, nbasis = P)
    psi.new[, l] <- psi.smooth$est.value
  }
  
  fname <- paste0('stan_fit/M2_whole_dp/dp1_', which(Tstart==starting.time), '.RData')
  load(fname)
  
  ## calculate conditional survival probabilities at T+\delta_T
  long.traj.1 <- long.traj.2 <- long.traj.3 <- long.traj.4 <- long.traj.5 <- 
    surv.traj <- data.frame(Time = t.new, mean = NA, q025 = NA, q975 = NA)
  long.traj <- list(long.traj.1, long.traj.2, long.traj.3, long.traj.4, long.traj.5)
  long.test.new <- long.test2[which(long.test2$ID==tmp.ID.new), ]
  surv.test.new <- surv.test2[which(surv.test2$ID==tmp.ID.new), ]
  
  time.index <- 0
  
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
    
    time.index <- time.index+1
    
    tmp.bs <- bs.mat.new[time.index, ]
    tmp.phi <- phi.new[time.index, ]
    tmp.psi <- psi.new[time.index, ]
    
    ## Obtain hazard index matrix ##
    t0 <- Tstart; t1 <- Tstart + Tdelta
    h.index.t0 <- h.index.t1 <- matrix(0, nrow = 1, ncol = n.tau)
    h.grid.t0 <- h.grid.t1 <- matrix(0, nrow = 1, ncol = n.tau)
    
    for (k in 1:n.tau){
      h.index.t0[1, k] <- ifelse(t0>=tau[k] & t0<tau[k+1], 1, 0)
      h.index.t1[1, k] <- ifelse(t1>=tau[k] & t1<tau[k+1], 1, 0)
      h.grid.t0[1, k] <- ifelse(t0>=tau[k], min(t0, tau[k+1]), t0)
      h.grid.t1[1, k] <- ifelse(t1>=tau[k], min(t1, tau[k+1]), t1)
    }
    h.grid.t0 <- cbind(0, h.grid.t0) 
    h.grid.t1 <- cbind(0, h.grid.t1)
    
    i <- which(surv.test2$ID==tmp.ID.new)
    xi_i <- xi_sample[, , i] 
    zeta1_i <- zeta1_sample[, , i]
    zeta2_i <- zeta2_sample[, , i]
    zeta3_i <- zeta3_sample[, , i]
    zeta4_i <- zeta4_sample[, , i]
    zeta5_i <- zeta5_sample[, , i]
    zeta_i <- list(zeta1_i, zeta2_i, zeta3_i, zeta4_i, zeta5_i)
    xi_m_i <- xi_m_sample[, i, ]
    
    for (j in 1:J){
      tmp.mu = sum(tmp.bs*A[[j]]) + beta[j]*(xi_i %*% tmp.phi + zeta_i[[j]]*tmp.psi)
      if (j==1) tmp.mu <- (tmp.mu*0.6+1)^(1/0.6)
      if (j==2) tmp.mu <- (tmp.mu*0.4+1)^(1/0.4)
      if (j==3) tmp.mu <- (tmp.mu*0.8+1)^(1/0.8)
      if (j==4) tmp.mu <- exp(tmp.mu)
      if (j==5) tmp.mu <- exp(tmp.mu)
      tmp.mu <- tmp.mu[which(!is.na(tmp.mu))]
      long.traj[[j]][time.index, ] <- c(t1, mean(tmp.mu), quantile(tmp.mu, c(0.025, 0.975)))
    }
    
    Hi_est_t0 <- Hi_est_t1 <- rep(NA, Q/2)
    for (q in (Q/2+1):Q){
      Hi_est_t0[q-Q/2] <- H_est(phi1_interval_T0, phi2_interval_T0, psi1_interval_T0, 
                                const_T0, x.test[i, ], weights, h_grid_T0, 
                                xi_i[q, ], zeta1_i[q], zeta2_i[q], zeta3_i[q], 
                                zeta4_i[q], zeta5_i[q], xi_m_i[q, ], 
                                t = Tstart)
      Hi_est_t1[q-Q/2] <- H_est(phi1_interval_T1, phi2_interval_T1, psi1_interval_T1, 
                                const_T1, x.test[i, ], weights, h_grid_T1, 
                                xi_i[q, ], zeta1_i[q], zeta2_i[q], zeta3_i[q], 
                                zeta4_i[q], zeta5_i[q], xi_m_i[q, ], 
                                t = Tstart + Tdelta)
    }
    cond_surv_prob <- exp(-Hi_est_t1 + Hi_est_t0)
    surv.traj[time.index, ] <- c(t1, mean(cond_surv_prob), 
                                 quantile(cond_surv_prob, c(0.025, 0.975)))
  }
  
  for (j in 1:J){
    long.traj[[j]] <- long.traj[[j]][which(long.traj[[j]]$Time<=surv.test.new$surv_time), ]
  }
  surv.traj <- surv.traj[which(surv.traj$Time<=surv.test.new$surv_time), ]
  
  ## 
  ylim <- list(c(0, 85), c(0, 75), c(-15, 15), c(0, 30), c(0, 18))
  xlim <- c(0, surv.test.new$surv_time)
  ylab <- c('ADAS-Cog 13', 'RAVLT-immediate', 'RAVLT-learning', 'MMSE', 'CDR-SB', 'Surv. Prob.')
  xlab <- c('', '', '', 'Time', 'Time', 'Time')
  long.test.new <- as.matrix(long.test.new)
  for (j in 1:J){
    y0 = long.test.new[, j+2]
    if (j==1) y0 <- (y0*0.6+1)^(1/0.6)
    if (j==2) y0 <- (y0*0.4+1)^(1/0.4)
    if (j==3) y0 <- (y0*0.8+1)^(1/0.8)
    if (j==4) y0 <- exp(y0)
    if (j==5) y0 <- exp(y0)
    plot(x = long.test.new[, 2], y = y0, type = 'p', 
         ylim = ylim[[j]], xlab = xlab[j], ylab = ylab[j], xlim = xlim)
    lines(x = long.traj[[j]]$Time, y = long.traj[[j]]$mean, lty = 1)
    lines(x = long.traj[[j]]$Time, y = long.traj[[j]]$q025, lty = 2)
    lines(x = long.traj[[j]]$Time, y = long.traj[[j]]$q975, lty = 2)
    abline(v = Tstart, lty = 2, col = 'blue')
  }
  
  plot(x = surv.traj$Time, y = surv.traj$mean, type = 'l', ylim = c(0, 1), 
       xlim = xlim, xlab = xlab[6], ylab = ylab[6])
  lines(x = surv.traj$Time, y = surv.traj$q025, lty = 2)
  lines(x = surv.traj$Time, y = surv.traj$q975, lty = 2)
  abline(v = Tstart, lty = 2, col = 'blue')
  abline(v = surv.test.new$surv_time, lty = 2-surv.test.new$status, col = 'red')
  
  dev.off()
}

}
