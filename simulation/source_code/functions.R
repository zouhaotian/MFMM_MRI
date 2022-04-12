align_mat <- function(Y, time, ID){
  uID <- length(unique(ID))
  argvals <- sort(unique(time))
  tg <- length(argvals)
  mat <- matrix(NA, uID, tg)
  
  for (tmp.ID in 1:uID){
    tmp.Y <- Y[which(ID==tmp.ID)]
    tmp.time <- time[which(ID==tmp.ID)]
    mat[tmp.ID, which(tmp.time %in% argvals)] <- tmp.Y
  }
  return(mat)
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


create_index <- function(time){
  S <- length(time)
  s.index <- rep(time, each = S)
  t.index <- rep(time, S)
  index.dat <- data.frame(s.index = s.index, t.index = t.index)
  return(index.dat)
}

bs.smooth <- function(value, argvals, argvals.new, nbasis = 9){
  bs.mat <- create_bs(argvals, argvals, nbasis = nbasis)
  bs.mat.pred <- create_bs(argvals, argvals.new, nbasis = nbasis)
  lm.obj <- lm(value ~ bs.mat - 1)
  coef.est <- unname(coef(lm.obj))
  est.value <- (bs.mat.pred %*% coef.est)[, 1]
  l <- list(coef.est = coef.est, est.value = est.value)
  return(l)
}

# create orthogonal B-spline matrix on the time grid #
create_bs <- function(time.grid, pred.time, nbasis = 9){
  knot.seq <- seq(min(time.grid), max(time.grid), length.out = nbasis - 2)
  knots <- expand.knots(knot.seq, order = 4)
  basis.fun <- SplineBasis(knots, order = 4)
  
  bs.mat <- matrix(NA, length(pred.time), nbasis)
  for (i in 1:length(pred.time)){
    bs.mat[i, ] <- evaluate(basis.fun, pred.time[i])
  }
  return(bs.mat)
}

estMean <- function(Y, time, time.grid, nbasis = 9){
  bs.mat <- create_bs(time.grid, time.grid, nbasis = nbasis)
  bs.mat.expand <- create_bs(time.grid, time, nbasis = nbasis)
  lm.obj <- lm(Y ~ bs.mat.expand - 1)
  beta_hat <- unname(lm.obj$coefficients)
  pred.values <- (bs.mat %*% beta_hat)[, 1]
  l <- list(mu_hat = pred.values, beta_hat = beta_hat)
  return(l)
}

estMean2 <- function(Y.mat, time.grid, nbasis = 9){
  mu_hat <- colMeans(Y.mat, na.rm = T)
  bs.mat <- create_bs(time.grid, time.grid, nbasis = nbasis)
  lm.obj <- lm(mu_hat ~ bs.mat - 1)
  beta_hat <- unname(lm.obj$coefficients)
  l <- list(mu_hat = mu_hat, beta_hat = beta_hat)
  return(l)
}

estBeta <- function(Y.list, ID, time, time.grid){
  J <- length(Y.list)
  tp <- length(time.grid)
  sigma_sq_hat <- rep(NA, J)
  beta_hat <- c(1, rep(NA, J-1))
  Cov_j_j <- vector('list', length = J)
  Cov_j_J <- vector('list', length = J-1)
  
  dat.mface <- list(y1 = data.frame(subj = ID, argvals = time, y = Y.list[[1]]),
                    y2 = data.frame(subj = ID, argvals = time, y = Y.list[[J]]))
  fit.mface <- mface.sparse(data = dat.mface, argvals.new = time.grid)
  Chat.new <- matrix(fit.mface$Chat.new@x, nrow = 2*tp)
  Chat.raw <- fit.mface$Chat.raw.diag.new
  
  C11_hat <- Chat.new[1:tp, 1:tp]
  CJJ_hat <- Chat.new[(tp+1):(2*tp), (tp+1):(2*tp)]
  C1J_hat <- Chat.new[1:tp, (tp+1):(2*tp)]
  C1J_hat <- (C1J_hat + t(C1J_hat))/2
  sigma_sq_hat[1] <- mean(Chat.raw[1:tp] - diag(C11_hat))
  sigma_sq_hat[J] <- mean(Chat.raw[(tp+1):(2*tp)] - diag(CJJ_hat))
  Cov_j_j[[1]] <- C11_hat
  Cov_j_j[[J]] <- CJJ_hat
  Cov_j_J[[1]] <- C1J_hat
  
  for (j in 2:(J-1)){
    dat.mface <- list(y1 = data.frame(subj = ID, argvals = time, y = Y.list[[j]]),
                      y2 = data.frame(subj = ID, argvals = time, y = Y.list[[J]]))
    fit.mface <- mface.sparse(data = dat.mface, argvals.new = time.grid)
    Chat.new <- matrix(fit.mface$Chat.new@x, nrow = 2*tp)
    Chat.raw <- fit.mface$Chat.raw.diag.new
    
    Cjj_hat <- Chat.new[1:tp, 1:tp]
    CjJ_hat <- Chat.new[1:tp, (tp+1):(2*tp)]
    CjJ_hat <- (CjJ_hat + t(CjJ_hat))/2
    beta_hat[j] <- reg(CjJ_hat, C1J_hat)
    sigma_sq_hat[j] <- mean(Chat.raw[1:tp] - diag(Cjj_hat))
    Cov_j_j[[j]] <- Cjj_hat
    Cov_j_J[[j]] <- CjJ_hat
  }
  
  dat.mface <- list(y1 = data.frame(subj = ID, argvals = time, y = Y.list[[J]]),
                    y2 = data.frame(subj = ID, argvals = time, y = Y.list[[2]]))
  fit.mface <- mface.sparse(data = dat.mface, argvals.new = time.grid)
  Chat.new <- matrix(fit.mface$Chat.new@x, nrow = 2*tp)
  Chat.raw <- fit.mface$Chat.raw.diag.new
  
  CJJ_hat <- Chat.new[1:tp, 1:tp]
  CJ2_hat <- Chat.new[1:tp, (tp+1):(2*tp)]
  CJ2_hat <- (CJ2_hat + t(CJ2_hat))/2
  
  dat.mface <- list(y1 = data.frame(subj = ID, argvals = time, y = Y.list[[1]]),
                    y2 = data.frame(subj = ID, argvals = time, y = Y.list[[2]]))
  fit.mface <- mface.sparse(data = dat.mface, argvals.new = time.grid)
  Chat.new <- matrix(fit.mface$Chat.new@x, nrow = 2*tp)
  Chat.raw <- fit.mface$Chat.raw.diag.new
  
  C12_hat <- Chat.new[1:tp, (tp+1):(2*tp)]
  C12_hat <- (C12_hat + t(C12_hat))/2
  beta_hat[J] <- reg(CJ2_hat, C12_hat)
  
  l <- list(sigma_sq_hat = sigma_sq_hat, beta_hat = beta_hat, 
            Cov_j_j = Cov_j_j, Cov_j_J = Cov_j_J)
  return(l)
}

estC0 <- function(Cov_j_J, beta_hat, time.grid, J){
  tp <- length(time.grid)
  C0_hat <- matrix(NA, tp, tp)
  
  for (s in 1:tp){
    for (t in 1:tp){
      tmp.Y <- tmp.X <- rep(NA, J-1)
      for (j in 1:(J-1)){
        tmp.Y[j] <- Cov_j_J[[j]][s, t]  
        tmp.X[j] <- beta_hat[j]*beta_hat[J]
      }
      lm.obj <- lm(tmp.Y ~ tmp.X - 1)
      C0_hat[s, t] <- unname(coef(lm.obj))
    }
  }
  return(C0_hat)
}

estC1 <- function(Cov_j_J, beta_hat, time.grid, J, C0_hat){
  tp <- length(time.grid)
  C1_hat <- matrix(NA, tp, tp)
  
  for (s in 1:tp){
    for (t in 1:tp){
      tmp.Y <- tmp.X <- rep(NA, J)
      for (j in 1:J){
        tmp.Y[j] <- Cov_j_j[[j]][s, t] - beta_hat[j]^2*C0_hat[s, t]
        tmp.X[j] <- beta_hat[j]^2
      }
      lm.obj <- lm(tmp.Y ~ tmp.X - 1)
      C1_hat[s, t] <- unname(coef(lm.obj))
    }
  }
  return(C1_hat)
}

## smoothed covariance matrix of Y: smooth diagonal using tensor-product spline
estCov <- function(Y.mat, mu_est){
  gp <- 1:ncol(Y.mat) ## grid points
  mu_expand <- matrix(rep(mu_est, nrow(Y.mat)), nrow = nrow(Y.mat), byrow = T)
  demeaned.dat <- Y.mat - mu_expand
  raw.cov.mat <- matrix(NA, nrow = ncol(Y.mat), ncol = ncol(Y.mat))
  
  for (i in gp){
    for (j in gp){
      tmp <- demeaned.dat[, i]*demeaned.dat[, j]
      raw.cov.mat[i, j] <- sum(tmp, na.rm = T)/(sum(!is.na(tmp)) - 1)
    }
  }
  
  s.index <- rep(gp, each = ncol(Y.mat))
  t.index <- rep(gp, ncol(Y.mat))
  
  cov2 <- raw.cov.mat
  diag(cov2) <- rep(NA, ncol(Y.mat))
  cov2.vec <- as.vector(cov2)
  index.dat <- data.frame(s.index = s.index, t.index = t.index)
  
  K.0 <- matrix(predict(gam(cov2.vec ~ te(s.index, t.index, k = 9), 
                            family = gaussian), index.dat), 
                ncol(Y.mat), ncol(Y.mat))
  K.0 <- (K.0 + t(K.0))/2
  sigma2_hat <- mean(diag(raw.cov.mat) - diag(K.0))
  
  l <- list(C = K.0, sigmasq_hat = sigma2_hat)
  return(l)
}

calc_cov <- function(Y1.mat, Y2.mat){
  gp <- 1:ncol(Y1.mat)
  raw.cov.mat <- matrix(NA, nrow = ncol(Y1.mat), ncol = ncol(Y1.mat))
  for (i in gp){
    for (j in gp){
      raw.cov.mat[i, j] <- cov(Y1.mat[, i], Y2.mat[, j], use = 'complete.obs')
    }
  }
  return(raw.cov.mat)
}

calc_cov2 <- function(mat){
  cov.mat <- matrix(0, nrow(mat), nrow(mat))
  for (i in 1:nrow(mat)){
    for (j in 1:nrow(mat)){
      cov.mat[i, j] <- cov(mat[i, ], mat[j, ])
    }
  }
  return(cov.mat)
}

reg <- function(C1, C2){
  C1.vec <- as.vector(C1)
  C2.vec <- as.vector(C2)
  lm.obj <- lm(C1.vec ~ C2.vec - 1)
  beta <- unname(lm.obj$coefficients)
  return(beta)
}


PCA <- function(C, time.grid, pve = 0.9, nbasis = 9, by = 1, L = NULL){
  eigen.obj <- eigen(C)
  values <- eigen.obj$values*by
  values.gt0 <- values[values>0]
  prop.var <- cumsum(values.gt0)/sum(values.gt0)
  if (is.null(L)){
    Lx <- min(which(prop.var>pve)) ## Explains >pve variance
  } else {
    Lx <- L
  }
  
  eigen.func <- eigen.obj$vectors[, 1:Lx]/sqrt(by)
  l <- list(values = values[1:Lx],
            eigen.func = eigen.func,
            pve = prop.var[Lx])
  return(l)
}

scale.eigen <- function(eigen.obj, L, scale.pars){
  values <- eigen.obj$values[1:L]
  eigen.func <- eigen.func.scaled <- eigen.obj$eigen.func[, 1:L]
  values.scaled <- values/scale.pars
  for (l in 1:L){
    eigen.func.scaled[, l] <- eigen.func[, l]*sqrt(scale.pars[l])
  }
  l <- list(values = values.scaled, 
            eigen.func = eigen.func.scaled)
}

sign_eigen <- function(eigen.func, index, sign.efunc){
  if (sign(eigen.func[index])==sign.efunc){
    return(eigen.func)
  } else {
    return(-eigen.func)
  }
}

calc_chol <- function(B, D){
  mat <- B %*% D %*% t(B)
  mat <- mat + diag(0.0001, nrow(mat))
  sigma <- sqrt(diag(mat))
  tmp.mat <- matrix(1, nrow(mat), ncol(mat))
  for (i in 1:nrow(mat)){
    for (j in 1:nrow(mat)){
      tmp.mat[i, j] <- mat[i, j]/(sigma[i]*sigma[j])
    }
  }
  chol.obj <- chol(tmp.mat)
  l <- list(sigma = sigma, L = t(chol.obj))
  return(l)
}

calcTol <- function(l1, l2){
  v1 <- unlist(l1)
  v2 <- unlist(l2)
  diff <- max(abs(v1 - v2))
  return(diff)
}

get_quantile_vector <- function(v){
  l <- rep(NA, 4)
  l[1] <- mean(v)
  l[2] <- sd(v)
  l[3] <- quantile(v, 0.025)
  l[4] <- quantile(v, 0.975)
  return(l)
}

get_quantile_matrix <- function(mat){
  tmp.mat <- matrix(NA, ncol(mat), 4)
  for (i in 1:ncol(mat)){
    tmp.mat[i, ] <- get_quantile_vector(mat[, i])
  }
  return(tmp.mat)
}

get_mean_function <- function(func.est){
  mean.func <- colMeans(func.est)
  return(mean.func)
}

matrix_to_vector <- function(func.est){
  tp <- nrow(func.est)
  L <- ncol(func.est)
  tmp.vec <- rep(NA, tp*L)
  for (i in 1:L){
    start <- (i-1)*tp + 1
    end <- i*tp
    tmp.vec[start:end] <- func.est[, i]
  }
  return(tmp.vec)
}

g.inits <- function(value){
  l <- length(value)
  ret.value <- rep(NA, l)
  for (i in 1:l){
    if (abs(value[i])<10){
      scale.value <- abs(value[i])/3
    } else {
      scale.value <- abs(value[i])/10
    }
    ret.value[i] <- value[i] + rnorm(1, 0, scale.value)
  }
  return(ret.value)
}


g.inits.sd <- function(value){
  l <- length(value)
  ret.value <- rep(NA, l)
  for (i in 1:l){
    ret.value[i] <- value[i] + runif(1, -value[i]/3, value[i]/3)
  }
  return(ret.value)
}

plot.f <- function(x.prior, y.prior, abline.axis, x.posterior, y.true, y.mean, y.lower, y.upper, x.limit, y.limit.min, y.limit.max, x.lab, y.lab){
  plot(x = x.prior, y = y.prior, xlim = x.limit, ylim = c(y.limit.min, y.limit.max), xlab = x.lab, ylab = y.lab)
  points(x = c(x.prior, x.posterior), y = c(y.prior, y.true), col = 'blue', pch = 19, cex = 1)
  lines(x = c(x.prior, x.posterior), y = c(y.prior, y.true), col = 'blue', lty = 1, lwd = 1)
  abline(v = abline.axis, lty = 2, lwd = 1)
  points(x = x.posterior, y = y.mean, col = 'red', pch = 2, cex = 1)
  lines(x = x.posterior, y = y.mean, col = 'red', lty = 2, lwd = 1)
  points(x = x.posterior, y = y.lower, col = 'red', pch = 1, cex = 1)
  lines(x = x.posterior, y = y.lower, col = 'red', lty = 2, lwd = 1)
  points(x = x.posterior, y = y.upper, col = 'red', pch = 1, cex = 1)
  lines(x = x.posterior, y = y.upper, col = 'red', lty = 2, lwd = 1)
}


## get estimated function and 95% CI, and plot #
plot.function <- function(time.grid, values, lab.x, lab.y, range.y){
  summary.f <- data.frame(time.grid = time.grid, values = values)
  p <- ggplot(data = summary.f, aes(x = time.grid, y = values)) + 
    geom_line(aes(y = values), linetype = 'solid', size = 0.7) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill=NA, size=1),
          plot.margin = margin(0.5, 0.25, 0, 0.25, "cm")) +
    xlab(TeX(lab.x)) +
    ylab(TeX(lab.y)) + 
    ylim(range.y[1], range.y[2])
  return(p)
}

# Approximate integral of f1(v)*f2(v) via trapezodial rule ##
c_int <- function(f1, f2, width){
  f3 = f1*f2
  v = length(f3)
  left_int = sum(f3[1:(v-1)]*width)
  right_int = sum(f3[2:v]*width)
  return((left_int+right_int)/2)
}

eigen2 <- function(mat, n.values = 10, V = 100, sign, index){
  res = eigs_sym(mat, k = n.values + 5, which = 'LM')
  l = list(values = res$values[1:n.values]/V, 
           vectors = res$vectors[, 1:n.values]*sqrt(V))
  for (i in 1:n.values){
    if (sign(l$vectors[index[i], i])!=sign[i]){
      l$vectors[, i] = -l$vectors[, i]
    }
  }
  return(l)
}

cov_approx <- function(cov.mat, V){
  x.index <- rep(1:V, V)
  y.index <- rep(1:V, rep(V, V))
  cov2 <- cov.mat
  diag(cov2) <- rep(NA, V)
  cov2.vec <- as.vector(cov2)
  index.dat <- data.frame(x = x.index, y = y.index, z = cov2.vec)
  xyz <- index.dat[complete.cases(index.dat), ]
  xy.est <- data.frame(x = x.index, y = y.index)
  mba.pts <- mba.points(xyz, xy.est)$xyz.est
  K.0 <- matrix(mba.pts[, 3], nrow = V, ncol = V)
  K.0 <- (K.0 + t(K.0))/2
  
  return(K.0)
}
