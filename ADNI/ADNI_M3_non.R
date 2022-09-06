## Computation of DIC, EAIC, EBIC, BF, CPO, LPML, LOOIC, WAIC for M3: non-mixed functional model ##
# setwd('/MMFPCA_MRI/ADNI/')
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
library(loo)
library(Amelia)
library(sn)
rstan_options(auto_write = FALSE)

set.seed(2020)

source('source_code/functions.R')
source('source_code/sFPCA.R')

long <- read.csv('dataset/ADNI_long_2.csv') ## _2
surv <- read.csv('dataset/ADNI_surv.csv')

N <- nrow(surv)
J <- 5 ## number of longitudinal outcomes
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
x <- as.matrix(cbind(surv$Age, surv$Gender, surv$Education, surv$APOE4))

p1 <- ggplot(data = baseline.cumulative.hazard, aes(x = time, y = hazard)) + 
  geom_line() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(face = "bold", color = "#993333"),
        axis.text.y = element_text(face = "bold", color = "#993333"),
        axis.line.x = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        axis.line.y = element_line(arrow = arrow(length = unit(0.3, "cm"))), 
        plot.title = element_text(hjust = 0.5)) + 
  xlab('Time (in years)') + 
  ylab('Estimated Cumulative Hazard Function') + 
  xlim(0, 8.5) + 
  ylim(0, 2.5) + 
  scale_x_continuous(breaks = c(0, 4, 6, 8)) + 
  ggtitle('Cumulative Hazard for ADNI')

save(list = 'p1', file = 'RData/Cumulative_Hazard_ADNI.RData')

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
md = stan_model('source_code/M3_non.stan')
stan_dat <- list(n = N, J = J, 
                 nobs = nrow(long), nmiss = Y_missing$len.missing, 
                 id_long = long$ID, 
                 L0 = L0, L1 = L1, P = P, P_surv = ncol(x), 
                 Y1 = Y.list.new[[1]], Y2 = Y.list.new[[2]], Y3 = Y.list.new[[3]],
                 Y4 = Y.list.new[[4]], Y5 = Y.list.new[[5]],
                 miss_index1 = miss.index[[1]], miss_index2 = miss.index[[2]], 
                 miss_index3 = miss.index[[3]], miss_index4 = miss.index[[4]], 
                 miss_index5 = miss.index[[5]],
                 time = long$Time, x = x, 
                 surv_time = surv$surv_time, status = surv$status,
                 b = bs.mat, phi = phi, psi = psi,
                 Ltau = n.tau, tau = tau[1:n.tau], 
                 h_grid = h.grid, h_index = h.index)

rnd <- 1/3
inits1 <- list(A1 = g.inits(mu_est[[1]]$coefficient), A2 = g.inits(mu_est[[2]]$coefficient),
               A3 = g.inits(mu_est[[3]]$coefficient), A4 = g.inits(mu_est[[4]]$coefficient), 
               A5 = g.inits(mu_est[[5]]$coefficient), 
               omega = g.inits.sd(sigma_hat),
               logh0 = g.inits(logh0_hat), 
               gamma_x = as.array(g.inits(gamma_x_hat)), 
               Y1_imp = mean_Y[[1]]+rnd, Y2_imp = mean_Y[[2]]+rnd, Y3_imp = mean_Y[[3]]+rnd, 
               Y4_imp = mean_Y[[4]]+rnd, Y5_imp = mean_Y[[5]]+rnd)
inits2 <- list(A1 = g.inits(mu_est[[1]]$coefficient), A2 = g.inits(mu_est[[2]]$coefficient),
               A3 = g.inits(mu_est[[3]]$coefficient), A4 = g.inits(mu_est[[4]]$coefficient), 
               A5 = g.inits(mu_est[[5]]$coefficient),
               omega = g.inits.sd(sigma_hat),
               logh0 = g.inits(logh0_hat), 
               gamma_x = as.array(g.inits(gamma_x_hat)), 
               Y1_imp = mean_Y[[1]]-rnd, Y2_imp = mean_Y[[2]]-rnd, Y3_imp = mean_Y[[3]]-rnd, 
               Y4_imp = mean_Y[[4]]-rnd, Y5_imp = mean_Y[[5]]-rnd)

inits <- list(c1 = inits1, c2 = inits2)
pars <- c('A1', 'A2', 'A3', 'A4', 'A5', 
          'omega', 
          'logh0', 'gamma_x',
          'LL', 
          'Y1_imp', 'Y2_imp', 'Y3_imp', 'Y4_imp', 'Y5_imp')
fitStan <- sampling(md, data = stan_dat, iter = 3000, warmup = 2000, 
                    chains = 2, thin=1, init = inits, pars = pars, seed = 2020,
                    control = list(adapt_delta = 0.8, max_treedepth=10))
fname <- paste0('stan_fit/', 'M3_non_posterior.RData')
save(list = 'fitStan', file = fname)
# load(fname)
summ <- summary(fitStan)$summary[1:(J*P + J + n.tau + ncol(x)), c(1, 3, 4, 8, 9, 10)]
write.csv(summ, file = 'RData/summary_posterior_M3_non.csv')

## Extract samples ##
Q <- 2000 ## samples
A_sample <- extract(fitStan, pars = c('A1', 'A2', 'A3', 'A4', 'A5'))
omega_sample <- extract(fitStan, par = 'omega')$omega
logh0_sample <- extract(fitStan, 'logh0')$logh0
gamma_x_sample <- extract(fitStan, pars = 'gamma_x')$gamma_x
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
    
    for (j in 1:J){
      mu_j <- sum(A_sample[[j]][q, ]*bs.mat[i, ])
      ll_long[q, i, j] <- dnorm(Y.list.tmp[[j]][i], mu_j, 
                                omega_sample[q, j],  log = T)
    }
    
    ll_long_full[q, ID[i]] <-  ll_long_full[q, ID[i]] + sum(ll_long[q, i, ])
  }
}

Dbar.long <- sum(-2*ll_long_full)/Q
Dbar.surv <- sum(-2*ll_surv)/Q
Dbar <- Dbar.long + Dbar.surv

np <- nrow(summ)
EAIC <- Dbar + 2*np
EBIC <- Dbar + log(N)*np

## Compute CPO, LPML ##
LPML.long <- 0
for (j in 1:J){
  ll_long_obs <- ll_long[, , j]
  CPO.long <- 1/apply(1/exp(ll_long_obs), 2, mean)
  LPML.long <- LPML.long + sum(log(CPO.long)) 
}
CPO.surv <- 1/apply(1/exp(ll_surv), 2, mean)
LPML.surv <- sum(log(CPO.surv))
LPML <- LPML.long + LPML.surv

## Compute Dhat ##
A <- lapply(1:J, FUN = function(x) colMeans(A_sample[[x]]))
omega <- apply(omega_sample, 2, 'mean')
logh0 <- apply(logh0_sample, 2, 'mean')
gamma_x <- apply(gamma_x_sample, 2, 'mean')
Y_imp <- lapply(1:J, FUN = function(x) colMeans(Y_imp_sample[[x]]))

Y.list.i.stan <- Y.list.o
for (j in 1:J){
  Y.list.i.stan[[j]][ miss.index[[j]] ] <- Y_imp[[j]]
}

## Calculation of log likelihood ##
ll_long_posterior <- 0
ll_surv_posterior <- 0
for (i in 1:nobs){
  for (j in 1:J){
    mu_j <- sum(A[[j]]*bs.mat[i, ])
    ll_long_posterior <- ll_long_posterior + dnorm(Y.list.i.stan[[j]][i], mu_j, 
                                                   omega[j], log = T)
  }
}

status <- surv$status
for (i in 1:N){
  
  logh0_obs = sum(h.index[i, ]*logh0)
  rs_w <- sum(gamma_x*x[i, ])
  h <- exp(logh0_obs + rs_w)
  H <- rep(NA, n.tau)
  for (k in 1:n.tau){
    H[k] = exp(rs_w)*exp(logh0[k])*(h.grid[i, k+1] - h.grid[i, k])
  }
  ll_surv_posterior <- ll_surv_posterior + status[i]*log(h) + sum(-H)
}

Dhat <- -2*(ll_long_posterior + ll_surv_posterior)
DIC <- 2*Dbar - Dhat

## Compute BF: P(y | M) ##
log_predictive_probability <- ll_long_posterior + ll_surv_posterior + (np/2)*log(2*pi)

## Compute LOOIC, WAIC ##
ll_full <- ll_long_full + ll_surv
rel_n_eff <- relative_eff(exp(ll_full), chain_id = rep(1:2, each = 1000))
looic <- loo(ll_full, r_eff = rel_n_eff, cores = 4)$estimates[3, 1]
waic <- waic(ll_full)$estimates[3, 1]

l <- list(Dbar.long = Dbar.long, Dbar.surv = Dbar.surv,
          np = np, EAIC = EAIC, EBIC = EBIC,
          LPML.long = LPML.long, LPML.surv = LPML.surv,LPML = LPML, 
          Dhat.long = ll_long_posterior, Dhat.surv = ll_surv_posterior,
          Dhat = Dhat, DIC = DIC,
          log_predictive_probability = log_predictive_probability,
          looic = looic, waic = waic)
save(list = 'l', file = 'RData/summary_posterior_M3_non_MP.RData')
