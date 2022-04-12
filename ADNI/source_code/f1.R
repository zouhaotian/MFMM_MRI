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


# Approximate integral of f1(v)*f2(v) via trapezodial rule ##
c_int <- function(f1, f2, width){
  f3 = f1*f2
  v = length(f3)
  left_int = sum(f3[1:(v-1)]*width)
  right_int = sum(f3[2:v]*width)
  return((left_int+right_int)/2)
}
