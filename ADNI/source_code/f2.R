expand_matrix <- function(ID, mat){
  count <- as.vector(table(ID))
  mat <- unname(mat)
  full.mat <- matrix(NA, nrow = length(ID), ncol = ncol(mat))
  full.mat <- full.mat[0, ]
  for (i in 1:nrow(mat)){
    tmp.mat <- matrix(rep(mat[i, ], count[i]), nrow = count[i], byrow = T)
    full.mat <- rbind(full.mat, tmp.mat)
  }
  return(full.mat)
}


expand_B <- function(dat, B){
  count <- as.vector(table(dat$ID))
  full.B <- matrix(0, nrow = sum(count), ncol = ncol(B))
  start <- 0; end <- 0
  for (i in 1:length(count)){
    start <- start + 1
    tmp.long <- dat[which(dat$ID==i), ]
    end <- end + count[i]
    tmp.B <- matrix(rep(B[i, ], count[i]), nrow = count[i], byrow = T)
    full.B[start:end, ] <- tmp.B
    start <- end
  }
  return(full.B)
}

calc_init_long1 <- function(ID, Y, time.spline, time, B){
  full.B <- expand_matrix(ID, B)
  lme.obj <- lmer(Y ~ time.spline + full.B + (1+time |ID))
  fixed <- lme.obj@beta
  re <- as.matrix(ranef(lme.obj)$ID)
  cov.est <- unname(VarCorr(lme.obj)$ID)
  var.est <- c(diag(cov.est), sigma(lme.obj)^2)
  corr.est <- cov.est[1, 2]/sqrt(cov.est[1, 1]*cov.est[2, 2])
  l <- list(fixed = fixed, re = re, var.est = var.est, corr.est = corr.est)
  return(l)
}

calc_init_long2 <- function(ID, Y, time.spline, re, B){
  full.B <- expand_matrix(ID, B)
  full.re <- expand_matrix(ID, re)
  lm.obj <- lm(Y ~ time.spline + full.B + full.re)
  fixed <- unname(lm.obj$coefficients)
  sigma_est <- sigma(lm.obj)^2
  l <- list(fixed = fixed, var.est = sigma_est)
  return(l)
}

index_determine <- function(ID){
  uID <- unique(ID)
  n <- length(uID)
  start <- end <- rep(0, n)
  for (i in 1:n){
    start[i] <- min(which(ID==uID[i]))
    end[i] <- max(which(ID==uID[i]))
  }
  l <- list(start = start, end = end)
  return(l)
}

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


create_B <- function(dat, mu_est, phi_est, obsgrid, L=8, by=0.01){
  Lx <- ncol(phi_est)
  N <- nrow(dat)
  S <- ncol(dat)
  bs_mat <- bs(obsgrid, df=L, intercept = T) ## create b-spline matrix
  
  int_mu_b <- mu_est %*% bs_mat * by
  int_phi_b <- t(phi_est) %*% bs_mat * by
  
  tmp.mu_est <- matrix(rep(as.vector(mu_est), N), nrow = N, byrow = T)
  ksi_est <- (dat - tmp.mu_est) %*% phi_est * by
  
  B.1 <- matrix(rep(as.vector(int_mu_b), N), byrow = T, nrow = N)
  B.2 <- ksi_est %*% int_phi_b
  B <- B.1 + B.2
  return(B)
}

