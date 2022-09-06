est_mean_cov <- function(){
  #### dataset information ####
  N = length(unique(long$ID))
  J = ncol(long) - 2 ## number of longitudinal outcomes
  ID = 1:N
  long.nobs = table(long$ID)
  ltnew = length(tnew)
  long.t = long[, 2]
  
  V <- ncol(MRI.dat)
  m.grid = (1:V)/V
  m.grid.full = rep(m.grid, N)
  
  #### Define true parameters ####
  if (!is.null(mu_m)){
    mu_m_true = array(NA, dim = c(V, ltm))
    for (v in 1:V){
      mu_m_true[v, ] = mu_m(v/V, tm)
    }
  } else mu_m_true = NULL
  
  if (is.null(phi_m.sign)){
    phi_m.sign = rep(1, L0)
    phi_m.index = rep(1, L0)
  }
  
  if (is.null(phi.sign)){
    phi.sign = rep(1, L0)
    phi.index = rep(1, L0)
  }
  
  ##### Estimate mu_m(v) and obtain demeaned MRI data ####
  tmp.MRI = NULL
  for (i in 1:N) tmp.MRI = c(tmp.MRI, MRI.dat[i, ])
  G = gam(tmp.MRI ~ s(m.grid.full, bs = 'ps'))
  MRI.v = data.frame(m.grid.full = m.grid)
  MRI.mu = predict(G, MRI.v)
  
  MRI.demean = MRI.dat
  for (i in 1:N){
    MRI.demean[i, ] <- MRI.dat[i, ] - MRI.mu
  }
  
  if (!is.null(mu_m_true)){
    MSE.mu_m = mean((MRI.mu - mu_m_true)^2)
  } else MSE.mu_m = NULL
  
  #### Estimate beta_m and phi_m(v) ####
  for (j in 1:J){
    Cmj_est = matrix(NA, nrow = ltnew, ncol = V)
    for (v in 1:V){
      tmp.MRI = MRI.demean[, v]
      tmp.MRI2 = rep(tmp.MRI, long.nobs)
      Cmj_raw = tmp.MRI2*long[, j+2]
      
      fit = gam(Cmj_raw ~ s(long.t, bs = 'ps'))
      new.dat2 = data.frame(long.t = tnew)
      Cmj_est[, v] <- predict(fit, new.dat2)
    }
    
    svd.obj <- svd(Cmj_est)
    tmp.phi <- svd.obj$u[, 1:L0]
    tmp.phi_m <- svd.obj$v[, 1:L0]
    
    for (l in 1:L0){
      tmp.phi_m[, l] <- sign_eigen(tmp.phi_m[, l], phi_m.index[l], phi_m.sign[l])
      tmp.phi_m[, l] <- tmp.phi_m[, l]*sqrt(V)
      tmp.phi[, l] <- sign_eigen(tmp.phi[, l], phi.index[l], phi.sign[l])
      tmp.phi[, l] <- tmp.phi[, l]*sqrt(ltnew)
    }
    
  }
  
  
  
  
  
}