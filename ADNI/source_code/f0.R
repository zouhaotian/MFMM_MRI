FPCA <- function(dat, obsgrid, pve=0.99, by=0.005, L=8){
  N <- nrow(dat)
  S <- ncol(dat)
  
  bs_mat <- bs(obsgrid, df = L, intercept = TRUE)
  X <- bs_mat
  Y <- colMeans(dat)
  lm.obj <- lm(Y ~ X - 1)
  coef <- lm.obj$coefficients %>% as.matrix()
  mu_est <- bs_mat %*% coef %>% as.vector()
  
  tmp.mu_est <- matrix(rep(mu_est, N), nrow = N, byrow = T)
  demean.dat <- dat - tmp.mu_est
  cov.mat.raw <- 1/(N-1) * t(demean.dat) %*% demean.dat ## raw covariance matrix ##
  s.index <- rep(0:(S-1), S)
  t.index <- rep(0:(S-1), rep(S, S))
  cov2 <- cov.mat.raw
  diag(cov2) <- rep(NA, S)
  cov2.vec <- as.vector(cov2)
  index.dat <- data.frame(s.index = s.index, t.index = t.index)
  
  K.0 <- matrix(predict(gam(as.vector(cov2.vec) ~ te(s.index, t.index, k=10)), index.dat), S, S)
  K.0 <- (K.0 + t(K.0))/2
  eigen.obj <- eigen(K.0)
  
  values <- eigen.obj$values*by
  values.gt0 <- values[values>0]
  prop.var <- cumsum(values.gt0)/sum(values.gt0)
  Lx <- min(which(prop.var>pve)) ## Explains >pve variance
  
  eigen.func <- eigen.obj$vectors[, 1:Lx]/sqrt(by)
  l <- list(mu_est = mu_est, phi_est = eigen.func, values = values.gt0[1:Lx],
            pve = prop.var[Lx])
  return(l)
}

FPCA_hd <- function(dat, obsgrid, pve = 0.9999, by=0.005, L=8){
  N <- nrow(dat)
  S <- ncol(dat)
  
  bs_mat <- bs(obsgrid, df = L, intercept = TRUE)
  X <- bs_mat
  Y <- colMeans(dat)
  lm.obj <- lm(Y ~ X - 1)
  coef <- lm.obj$coefficients %>% as.matrix()
  mu_est <- bs_mat %*% coef %>% as.vector()
  
  tmp.mu_est <- matrix(rep(mu_est, N), nrow = N, byrow = T)
  demean.dat <- dat - tmp.mu_est
  cov.mat.raw <- 1/(N-1) * t(demean.dat) %*% demean.dat ## raw covariance matrix ##
  x.index <- rep(0:(S-1), S)
  y.index <- rep(0:(S-1), rep(S, S))
  cov2 <- cov.mat.raw
  diag(cov2) <- rep(NA, S)
  cov2.vec <- as.vector(cov2)
  index.dat <- data.frame(x = x.index, y = y.index, z = cov2.vec)
  xyz <- index.dat[complete.cases(index.dat), ]
  xy.est <- data.frame(x = x.index, y = y.index)
  mba.pts <- mba.points(xyz, xy.est)$xyz.est
  K.0 <- matrix(mba.pts[, 3], nrow = S, ncol = S)
  K.0 <- (K.0 + t(K.0))/2
  
  res = eigs_sym(K.0, k = 100, which = "LM")
  
  values <- res$values*by
  values.gt0 <- values[values>0]
  prop.var <- cumsum(values.gt0)/sum(values.gt0)
  Lx <- min(which(prop.var>pve)) ## Explains >pve variance
  
  eigen.func <- res$vectors[, 1:Lx]/sqrt(by)
  l <- list(mu_est = mu_est, phi_est = eigen.func, values = values.gt0[1:Lx],
            pve = prop.var[Lx])
  return(l)
}

