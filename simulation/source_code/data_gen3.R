## data simulation function ##
## N: number of subjects; seed: random seed ##
sim <- function(seed, tp, N){
  set.seed(seed)
  
  J <- tp$J
  V <- tp$V
  m.grid <- (1:V)/V
  tnew <- tp$tnew
  
  mu <- tp$mu
  
  mu_m <- function(v) 1 + sin(2*pi*v) + cos(2*pi*v)
  mu_m_obs = mu_m(m.grid)
  
  L0 <- tp$L0
  phi <- list(
    f1 = function(t) sqrt(2)*sin(pi*t),
    f2 = function(t) sqrt(2)*cos(pi*t)
  )
  
  L1 <- tp$L1
  psi <- list(
    f1 = function(t) sqrt(2)*cos(2*pi*t)
  )
  
  phi_obs <- matrix(NA, nrow = length(tnew), ncol = L0)
  psi_obs <- matrix(NA, nrow = length(tnew), ncol = L1)
  for (l in 1:L0) phi_obs[, l] <- phi[[l]](tnew)
  for (l in 1:L1) psi_obs[, l] <- psi[[l]](tnew)
  
  ## Eigenfunction of u_mi(v)
  phi_m <- list(
    f1 = function(v) sqrt(2)*sin(pi*v/2),
    f2 = function(v) sqrt(2)*sin(3*pi*v/2)
  )
  
  ## Eigenfunction of f_mi(v)
  Lm <- tp$Lm
  psi_m <- list(
    f1 = function(v) sqrt(2)*cos(pi*v/2),
    f2 = function(v) sqrt(2)*cos(3*pi*v/2),
    f3 = function(v) sqrt(2)*cos(5*pi*v/2),
    f4 = function(v) sqrt(2)*cos(7*pi*v/2),
    f5 = function(v) sqrt(2)*cos(9*pi*v/2),
    f6 = function(v) sqrt(2)*cos(11*pi*v/2)
  )
  
  phi_m_obs <- matrix(NA, nrow = length(m.grid), ncol = L0)
  psi_m_obs <- matrix(NA, nrow = length(m.grid), ncol = Lm)
  for (l in 1:L0) phi_m_obs[, l] <- phi_m[[l]](m.grid)
  for (l in 1:Lm) psi_m_obs[, l] <- psi_m[[l]](m.grid)
  
  d0 <- tp$d0
  d1 <- tp$d1
  d_m <- tp$d_m
  
  beta <- tp$beta
  beta_m <- tp$beta_m
  omega <- tp$omega
  omega_m <- tp$omega_m
  
  C0 = phi_obs %*% diag(d0) %*% t(phi_obs)
  C1 = psi_obs %*% diag(d1, nrow = L1) %*% t(psi_obs)
  Ct = C0 + C1
  
  K.1 = beta_m^2*(phi_m_obs %*% diag(d0) %*% t(phi_m_obs))
  K.2 = psi_m_obs %*% diag(d_m) %*% t(psi_m_obs)
  cov_mat_m = K.1 + K.2
  
  logh0 <- tp$logh0
  gamma_x <- tp$gamma_x
  gamma0 <- c(tp$gamma01, tp$gamma02)
  gamma11 <- tp$gamma11
  gamma12 <- tp$gamma12
  gamma13 <- tp$gamma13
  gamma_m <- tp$gamma_m
  
  xi <- matrix(NA, N, L0)
  zeta1 <- zeta2 <- zeta3 <- matrix(NA, N, L1)
  xi_m <- matrix(NA, N, Lm)
  
  for (l in 1:L0){
    xi[, l] <- (sqrt(d0[l])*scale(rnorm(N)))[, 1]
  }
  
  for (l in 1:L1){
    zeta1[, l] <- (sqrt(d1[l])*scale(rnorm(N)))[, 1]
    zeta2[, l] <- (sqrt(d1[l])*scale(rnorm(N)))[, 1]
    zeta3[, l] <- (sqrt(d1[l])*scale(rnorm(N)))[, 1]
  }
  zeta <- list(zeta1, zeta2, zeta3)
  
  for (l in 1:Lm){
    xi_m[, l] <- (sqrt(d_m[l])*scale(rnorm(N)))[, 1]
  }
  
  int_m_psi <- matrix(NA, nrow = N, ncol = Lm)
  for (i in 1:N){
    mi = as.function(alist(v = , beta_m*(xi[i, 1]*phi_m[[1]](v) + xi[i, 2]*phi_m[[2]](v)) + 
                             xi_m[i, 1]*psi_m[[1]](v) + xi_m[i, 2]*psi_m[[2]](v) + xi_m[i, 3]*psi_m[[3]](v) + 
                             xi_m[i, 4]*psi_m[[4]](v) + xi_m[i, 5]*psi_m[[5]](v) ))
    for (j in 1:Lm){
      f <- as.function(alist(v = , mi(v)*psi_m[[j]](v)))
      int_m_psi[i, j] <- stats::integrate(f, lower=0, upper=1)$value
    }
  }
  
  ID <- time <- time.points <- NULL
  mu_obs <- vector('list', J)
  phi_obs <- vector('list', L0)
  psi_obs <- vector('list', L1)
  for (i in 1:N){
    tmp.obstime.2 <- sample((1:10)/100, size = 1)
    tmp.obstime.2 <- tmp.obstime.2 + (0:10)/10
    if (tmp.obstime.2[11]>max(tnew)) tmp.obstime.2 <- tmp.obstime.2[-11]
    tmp.obstime.2 <- c(0, tmp.obstime.2)
    
    time <- c(time, tmp.obstime.2)
    time.points <- c(time.points, length(tmp.obstime.2))
    ID <- c(ID, rep(i, length(tmp.obstime.2)))
    time.index = tmp.obstime.2*100 + 1
    
    for (j in 1:J) mu_obs[[j]] <- c(mu_obs[[j]], mu[[j]][time.index])
    
    for (l in 1:L0) phi_obs[[l]] <- c(phi_obs[[l]], phi[[l]](tmp.obstime.2))
    for (l in 1:L1) psi_obs[[l]] <- c(psi_obs[[l]], psi[[l]](tmp.obstime.2))
  }
  
  s.time.points <- sum(time.points)
  phi_obs_mat <- matrix(NA, s.time.points, L0)
  psi_obs_mat <- matrix(NA, s.time.points, L1)
  for (l in 1:L0) phi_obs_mat[, l] <- phi_obs[[l]]
  for (l in 1:L1) psi_obs_mat[, l] <- psi_obs[[l]]
  
  Y <- true.Y <- err <- vector('list', J)
  for (j in 1:J) err[[j]] <- (omega[j]*scale(rnorm(s.time.points)))[, 1]
  
  for (i in 1:s.time.points){
    for (j in 1:J){
      tmp.Y <- mu_obs[[j]][i] + beta[j]*(sum(xi[ID[i], ]*phi_obs_mat[i, ]) + 
                                           sum(zeta[[j]][ID[i], ]*psi_obs_mat[i, ]))
      true.Y[[j]] <- c(true.Y[[j]], tmp.Y)
      Y[[j]] <- c(Y[[j]], tmp.Y + err[[j]][i])
    }
  }
  
  long.dat <- data.frame(ID = ID, Time = time, 
                         y1 = Y[[1]], y2 = Y[[2]], y3 = Y[[3]])
  
  m_mat <- err_mat <- m_true <- matrix(NA, nrow = N, ncol = V)
  u_m <- f_m <- matrix(NA, nrow = N, ncol = V)
  for (i in 1:N){
    u_mi = (xi[i, ] %*% t(phi_m_obs))[1, ]
    f_mi = (xi_m[i, ] %*% t(psi_m_obs))[1, ]
    err_mi <- (omega_m*scale(rnorm(V)))[, 1]
    err_mat[i, ] <- err_mi
    m_mat[i, ] <- mu_m_obs + beta_m*u_mi + f_mi + err_mi
    m_true[i, ] <- mu_m_obs + beta_m*u_mi + f_mi
    u_m[i, ] <- u_mi
    f_m[i, ] <- f_mi
  }
  
  C <- runif(N, 0, 1.6)   ## censoring time
  C <- pmin(C, rep(max(tnew), N))
  surv.dat <- data.frame(ID = rep(1:N), x = rbinom(N, 2, 0.4), 
                         surv_time = rep(NA, N),
                         status = rep(NA, N))
  
  for (i in 1:N){
    Si <- runif(1, 0, 1) ## survival probability
    rs_w = sum(surv.dat$x[i]*gamma_x) + sum(gamma0*xi[i, ]) + 
      sum(gamma11*zeta1[i, ]) + sum(gamma12*zeta2[i, ]) + sum(gamma13*zeta3[i, ]) + 
      sum(gamma_m*xi_m[i, ])
    ti = -log(Si)/exp(logh0 + rs_w)
    surv.dat$surv_time[i] <- min(ti, C[i]) ## observed survival time
    surv.dat$status[i] <- as.integer(ti<C[i])
  }
  
  long.dat2 <- long.dat[0, ]
  for (i in 1:N){
    tmp.long.dat <- long.dat[which(long.dat$ID==i), ]
    tmp.surv.time <- surv.dat$surv_time[i]
    index <- which(tmp.long.dat$Time<=tmp.surv.time)
    long.dat2 <- rbind(long.dat2, tmp.long.dat[index, ])
  }
  
  ## Generate missingness ##
  inv_logit <- function(x) return(1/(1+exp(-x)))
  long.dat3 <- long.dat2
  for (i in 1:N){
    p_missing <- inv_logit(0.5*(long.dat2$Time[i]-5))
    I_missing <- rbinom(3, 1, p_missing)
    if (I_missing[1]==1) long.dat3$y1[i] <- NA
    if (I_missing[2]==1) long.dat3$y2[i] <- NA
    if (I_missing[3]==1) long.dat3$y3[i] <- NA
  }
  
  sim.dat <- list(long = long.dat3, surv = surv.dat, 
                  MRI.dat = m_mat, 
                  xi = xi, zeta = zeta, xi_m = xi_m, 
                  err_mat = err_mat, int_m_psi = int_m_psi, 
                  cov_mat_m = cov_mat_m, K.1 = K.1, K.2 = K.2, 
                  M_true = m_true, u_m = u_m, f_m = f_m, 
                  C0 = C0, C1 = C1, Ct = Ct)
  return(sim.dat)
}
