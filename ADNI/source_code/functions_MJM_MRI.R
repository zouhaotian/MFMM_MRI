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
  lme.obj <- lmer(Y ~ time.spline + full.B +  (1+time |ID))
  fixed <- lme.obj@beta
  re <- as.matrix(ranef(lme.obj)$ID)
  cov.est <- unname(VarCorr(lme.obj)$ID)
  var.est <- c(diag(cov.est), sigma(lme.obj)^2)
  corr.est <- cov.est[1, 2]/sqrt(cov.est[1, 1]*cov.est[2, 2])
  l <- list(fixed = fixed, re = re, var.est = var.est, corr.est = corr.est)
  return(l)
}

calc_init_long2 <- function(ID, Y, time.spline, B, re){
  full.B <- expand_matrix(ID, B)
  full.re <- expand_matrix(ID, re)
  lm.obj <- lm(Y ~ time.spline + full.B + full.re)
  fixed <- unname(lm.obj$coefficients)
  sigma_est <- sigma(lm.obj)^2
  l <- list(fixed = fixed, var.est = sigma_est)
  return(l)
}

calc_init_long1_without_B <- function(ID, Y, time.spline, time){
  lme.obj <- lmer(Y ~ time.spline + (1+time |ID))
  fixed <- lme.obj@beta
  re <- as.matrix(ranef(lme.obj)$ID)
  cov.est <- unname(VarCorr(lme.obj)$ID)
  var.est <- c(diag(cov.est), sigma(lme.obj)^2)
  corr.est <- cov.est[1, 2]/sqrt(cov.est[1, 1]*cov.est[2, 2])
  l <- list(fixed = fixed, re = re, var.est = var.est, corr.est = corr.est)
  return(l)
}

calc_init_long2_without_B <- function(ID, Y,  time.spline, re){
  full.re <- expand_matrix(ID, re)
  lm.obj <- lm(Y ~ time.spline + full.re)
  fixed <- unname(lm.obj$coefficients)
  sigma_est <- sigma(lm.obj)^2
  l <- list(fixed = fixed, var.est = sigma_est)
  return(l)
}


calc_init_surv <- function(surv_time, status, W.mat, B){
  cox.obj <- coxph(Surv(surv_time, status) ~ W.mat + B, method = 'breslow')
  paras.est <- cox.obj$coefficients %>% as.vector()
  return(paras.est)
}

calc_init_surv_withoutB <- function(surv_time, status, W.mat){
  cox.obj <- coxph(Surv(surv_time, status) ~ W.mat, method = 'breslow')
  paras.est <- cox.obj$coefficients %>% as.vector()
  return(paras.est)
}

missing_determine <- function(value){
  len <- length(value)
  new.value <- vector('list', length = len)
  missing.index <- vector('list', length = len)
  len.missing <- rep(NA, len)
  for (j in 1:len){
    tmp.value <- value[[j]]
    tmp.missing.index <- which(is.na(tmp.value))
    tmp.value[tmp.missing.index] <- -100
    new.value[[j]] <- tmp.value
    missing.index[[j]] <- tmp.missing.index
    len.missing[j] <- length(tmp.missing.index)
  }
  l <- list(new.value = new.value, missing.index = missing.index,
            len.missing = len.missing)
}


missing_determine_new <- function(value, value.imp){
  len <- length(value)
  new.value <- vector('list', length = len)
  missing.index <- vector('list', length = len)
  len.missing <- rep(NA, len)
  for (j in 1:len){
    tmp.value <- value[[j]]
    tmp.missing.index <- which(is.na(tmp.value))
    tmp.value[tmp.missing.index] <- value.imp[[j]]
    new.value[[j]] <- tmp.value
    missing.index[[j]] <- tmp.missing.index
    len.missing[j] <- length(tmp.missing.index)
  }
  l <- list(new.value = new.value, missing.index = missing.index,
            len.missing = len.missing)
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

filter.time <- function(long, surv, t.start){
  surv.filtered <- surv[surv$surv_time>t.start, ]
  long.filtered <- long[long$ID %in% surv.filtered$ID, ]
  l <- list(surv.filtered = surv.filtered, long.filtered = long.filtered)
  return(l)
}

get_h <- function(curr.time){
  tmp.h.grid <- rep(0, n.tau)
  for (k in 1:n.tau){
    tmp.h.grid[k] <- ifelse(curr.time>=tau[k], min(curr.time, tau[k+1]), curr.time)
  }
  tmp.h.grid <- c(0, tmp.h.grid)
  return(tmp.h.grid)
}

get_quantile <- function(x){
  v <- rep(0, 3)
  v[1] <- mean(x)
  v[2] <- quantile(x, 0.025)
  v[3] <- quantile(x, 0.975)
  return(v)
}

plot1 <- function(ft, fY, ot, oY, rt, rY, Tstart, main, ylim, xlab, ylab){
  plot(ft, fY[, 1], type = 'l', col = 'blue', xlim = c(0, 14), ylim = ylim, main = main, xlab = xlab, ylab = ylab)
  lines(ft, fY[, 2], lty = 2)
  lines(ft, fY[, 3], lty = 2)
  abline(v = Tstart, lty = 3)
  points(ot, oY, pty = 2)
  points(rt, rY, pty = 1)
}

plot2 <- function(ft, fY, rt, rstatus, Tstart, main, ylim, xlab, ylab){
  plot(ft, fY[, 1], type = 'l', col = 'blue', xlim = c(0, 14), ylim = ylim, main = main, xlab = xlab, ylab = ylab)
  lines(ft, fY[, 2], lty = 2)
  lines(ft, fY[, 3], lty = 2)
  abline(v = Tstart, lty = 3)
  abline(v = rt, lty = 2 - rstatus, col = 'red')
}

get_tau <- function(curr.time){
  tmp.tau <- rep(0, n.tau)
  for (k in 1:n.tau){
    tmp.tau[k] <- ifelse(curr.time>=tau[k], curr.time - tau[k], 0)
  }
  return(tmp.tau)
}

get_tau_index <- function(curr.time){
  tmp.tau.index <- rep(0, n.tau)
  for (k in 1:n.tau){
    tmp.tau.index[k] <- ifelse(curr.time>=tau[k], 1, 0)
  }
  return(tmp.tau.index)
}

get_tau_const <- function(curr.time){
  tmp.tau.const <- rep(0, n.tau)
  for (k in 1:n.tau){
    tmp.tau.const[k] <- ifelse(curr.time>=tau[k], -tau[k], 0)
  }
  return(tmp.tau.const)
}

get_list <- function(l, ID, train = T){
  tmp.ID <- l$ID
  tmp.ID.new <- tmp.ID[which(tmp.ID %in% ID)]
  tmp.time <- l$time[which(tmp.ID %in% ID)]
  tmp.value <- l$value[which(tmp.ID %in% ID)]
  
  if (train==T) tmp.ID.new <- rep(1:length(ID), table(tmp.ID.new))
  l.train <- list(ID = tmp.ID.new, time = tmp.time, value = tmp.value)
  return(l.train)
}

