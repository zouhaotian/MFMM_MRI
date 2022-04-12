## M1: training MRI dataset; M2: testing MRI dataset
## xi1: estimated xi for training dataset; xi2: estimated xi for testing dataset
## L0: FPCs for u_{mi}(v); Lm: FPCs for f_mi(v); V: number of voxels
## phi_m.sign: sign of eigenfunction \phi_ml(v) for the index of \phi_ml(v)
## psi_m.sign: sign of eigenfunction \psi_ml(v) for the index of \psi_ml(v)
## max.iters: Maximum number of iterations; 
## eps: maximum epsilon allowed (difference in \beta_m between new values and old values)
## beta_m, Phi_m_est: initial values of beta_m and Phi_m_est

search <- function(M1, M2, xi1, xi2, L0, Lm, V, 
                   phi_m.sign, phi_m.index, 
                   psi_m.sign, psi_m.index, 
                   max.iters = 20, eps = 1e-5, 
                   beta_m = NULL, Phi_m_est = NULL, 
                   d0, sign_beta_m){
  
  iters = 0    ## current iteration
  curr.eps = 100  ## current epsilon
  curr.eps.phi <- NULL
  beta_m_all <- NULL
  
  cov.M1 = cov(M1)
  
  ## Initial values of beta_m and Phi_all
  if (is.null(beta_m) | is.null(Phi_m_est)){
    cov.mat.approx = cov_approx(cov.M1, V)
    res = eigen2(cov.mat.approx, n.values = L0, V = V, sign = phi_m.sign, index = phi_m.index)
    
    X = d0
    Y = res$values
    beta_m = sign_beta_m*sqrt(1/sum(X*X)*sum(X*Y))
    Phi_m_est = res$vector
  }
  
  while (iters<max.iters & curr.eps>eps){
    beta_m_all = c(beta_m_all, beta_m)
    iters = iters + 1
    
    ## Estimate \psi_{ml}(v) and \xi_{mil} ##
    xi_phi_m = beta_m*(xi2 %*% t(Phi_m_est))
    tmp.M = M2 - xi_phi_m
    cov.mat = cov(tmp.M)
    cov.mat.approx = cov_approx(cov.mat, V)
    
    res = eigen2(cov.mat.approx, n.values = Lm, V = V, sign = psi_m.sign, index = psi_m.index)
    Psi_m_est = res$vector
    dm = res$values
    
    ## Estimate \phi_{ml}(v) and \beta_m ##
    cov.mat = cov.M1 - Psi_m_est %*% diag(dm) %*% t(Psi_m_est)
    cov.mat.approx = cov_approx(cov.mat, V)
    res = eigen2(cov.mat.approx, n.values = L0, V = V, sign = phi_m.sign, index = phi_m.index)
    Phi_m_est_n = res$vector
    
    X = d0
    Y = res$values
    beta_m_n = sign_beta_m*sqrt(1/sum(X*X)*sum(X*Y))
    
    curr.eps = abs(beta_m_n - beta_m)
    beta_m = beta_m_n
    
    tmp.eps.phi = mean(abs(Phi_m_est_n - Phi_m_est)^2)
    curr.eps.phi = c(curr.eps.phi, tmp.eps.phi)
    Phi_m_est = Phi_m_est_n
  }
  
  if (iters>=max.iters){
    warning("Reached maximum number of iterations without convergence!")
  } else {
    cat(sprintf("After %d iterations, model converged!\n", iters))
  }
  
  l = list(beta_m = beta_m, Phi_m_est = Phi_m_est, Psi_m_est = Psi_m_est, 
           dm = dm, 
           beta_m_all = beta_m_all, curr.eps.phi = curr.eps.phi)
}

plot.function <- function(obsgrid, f.index, f.true, f.all, ylim){
  plot(obsgrid, f.true[, f.index], type = 'n', ylim = ylim)
  lines(obsgrid, f.true[, f.index], lwd = 2, lty = 1)
  
  for (i in 1:length(f.all)){
    lines(obsgrid, f.all[[i]][, f.index], lwd = 1, lty = 2, col = 'red')
  }
}
