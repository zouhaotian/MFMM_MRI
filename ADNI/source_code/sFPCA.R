
## post-process covariance
sFPCA_fit <- function(S_fit, pve = 0.99){
  eigS = eigen(S_fit)
  sel = (eigS$values > 0)
  eigvalues = eigS$values[sel]
  eigvectors = eigS$vectors[,sel]
  npc = which.max(cumsum(eigvalues) / sum(eigvalues) > pve)
  eigvalues = eigvalues[1:npc] # <---
  eigvectors = eigvectors[,1:npc] # <---
  if (npc==1){
    S_fit = eigvectors %*% t(eigvectors) * eigvalues
  }else{
    S_fit = eigvectors %*% diag(eigvalues) %*% t(eigvectors)
  }
  return(list(S = S_fit, eigvalues = eigvalues, eigvectors = eigvectors, 
              pve = pve, npc = npc))
}


sFPCA_post_all <- function(S0_fit, S1_fit, 
                           mfaces_fit, Beta, pve = 0.99){
  
  p <- 5
  
  S0 = S0_fit$S
  S1 = S1_fit$S
  
  beta1 = Beta[1]
  beta2 = Beta[2]
  beta3 = Beta[3]
  beta4 = Beta[4]
  beta5 = Beta[5]
  
  fit1 = mfaces_fit$fit$y1
  fit2 = mfaces_fit$fit$y2
  fit3 = mfaces_fit$fit$y3
  fit4 = mfaces_fit$fit$y4
  fit5 = mfaces_fit$fit$y5
  
  object = mfaces_fit
  tnew = object$argvals.new
  Bnew <- vector("list", p)
  B <- vector("list", p)
  for(k in 1:p){
    Bnew[[k]] <- object$fit[[k]]$Bnew %*% object$G_invhalf
    B[[k]] <- object$fit[[k]]$B %*% object$G_invhalf
  }
  Bbasis = object$fit[[k]]$Bnew %*% object$G_invhalf
  Bn = Bnew[[1]] / nrow(object$var.error.new)
  
  Theta0 = t(Bn) %*% S0 %*% Bn
  Theta1 = t(Bn) %*% S1 %*% Bn
  # Theta2 = t(Bn) %*% S2 %*% Bn
  # Theta3 = t(Bn) %*% S3 %*% Bn
  # Theta4 = t(Bn) %*% S4 %*% Bn
  # Theta5 = t(Bn) %*% S5 %*% Bn
  
  Theta11 = beta1^2*(Theta0 + Theta1)
  Theta12 = beta1*beta2*Theta0
  Theta13 = beta1*beta3*Theta0
  Theta14 = beta1*beta4*Theta0
  Theta15 = beta1*beta5*Theta0
  
  Theta21 = beta2*beta1*Theta0
  Theta22 = beta2^2*(Theta0 + Theta1)
  Theta23 = beta2*beta3*Theta0
  Theta24 = beta2*beta4*Theta0
  Theta25 = beta2*beta5*Theta0
  
  Theta31 = beta3*beta1*Theta0
  Theta32 = beta3*beta2*Theta0
  Theta33 = beta3^2*(Theta0 + Theta1)
  Theta34 = beta3*beta4*Theta0
  Theta35 = beta3*beta5*Theta0
  
  Theta41 = beta4*beta1*Theta0
  Theta42 = beta4*beta2*Theta0
  Theta43 = beta4*beta3*Theta0
  Theta44 = beta4^2*(Theta0 + Theta1)
  Theta45 = beta4*beta5*Theta0
  
  Theta51 = beta5*beta1*Theta0
  Theta52 = beta5*beta2*Theta0
  Theta53 = beta5*beta3*Theta0
  Theta54 = beta5*beta4*Theta0
  Theta55 = beta5^2*(Theta0 + Theta1)
  
  
  Theta_all = rbind(cbind(Theta11, Theta12, Theta13, Theta14, Theta15),
                    cbind(Theta21, Theta22, Theta23, Theta24, Theta25),
                    cbind(Theta31, Theta32, Theta33, Theta34, Theta35),
                    cbind(Theta41, Theta42, Theta43, Theta44, Theta45),
                    cbind(Theta51, Theta52, Theta53, Theta54, Theta55))
  Bnew <- do.call(bdiag, Bnew)#%*%G_invhalfmat
  Chat <- Bnew%*%Theta_all%*%t(Bnew)
  
  Theta_W = bdiag(beta1^2*Theta1, beta2^2*Theta1, beta3^2*Theta1, beta4^2*Theta1, beta5^2*Theta1)
  Theta_U = Theta_all - Theta_W
  
  
  Chat_diag = as.vector(diag(Chat))  
  Cor = diag(1/sqrt(Chat_diag))%*%Chat%*%diag(1/sqrt(Chat_diag))
  
  B <- do.call(bdiag, B)
  Chat.diag.pred = diag(as.matrix(tcrossprod(B%*%Matrix(Theta_all),B)))
  
  Eig <- eigen(Theta_all)
  npc <- which.max(cumsum(Eig$values)/sum(Eig$values)>pve)[1]
  
  eigenfunctions = matrix(Bnew%*%Eig$vectors[,1:min(npc, p*length(tnew))],
                          ncol=min(npc, p*length(tnew)))
  eigenvalues = Eig$values[1:min(npc, p*length(tnew))]
  
  Eig0 <- eigen(Theta0)
  npc0 <- which.max(cumsum(Eig0$values)/sum(Eig0$values)>pve)[1]
  
  eigenfunctions0 = matrix(Bbasis%*%Eig0$vectors[,1:min(npc0, length(tnew))],
                           ncol=min(npc0, length(tnew)))
  eigenvalues0 = Eig0$values[1:min(npc0, length(tnew))]
  
  Eig1 <- eigen(Theta1)
  npc1 <- which.max(cumsum(Eig1$values)/sum(Eig1$values)>pve)[1]
  
  eigenfunctions1 = matrix(Bbasis%*%Eig1$vectors[,1:min(npc1, length(tnew))],
                           ncol=min(npc1, length(tnew)))
  eigenvalues1 = Eig1$values[1:min(npc1, length(tnew))]
  
  # r_j = matrix(NaN, nrow(data),p)
  # for(j in 1:p){
  #   obs_j = data[,-c(1,2)][,j]
  #   r_j[,j] = (obs_j - object$fit[[j]]$mu.hat)^2
  # }
  # 
  # sigma2 = colMeans(matrix(c(r_j) - Chat.diag.pred, nrow(data), p))
  # sigma2[sigma2 <= 0] <- 1e-6
  # 
  # var.error.hat <- t(matrix(rep(sigma2, p*nrow(data)),p,nrow(data)))
  # var.error.new <- t(matrix(rep(sigma2, p*length(tnew)),p,length(tnew)))
  # 
  # for(j in 1:p){
  #   object$fit[[j]]$sigma2 = sigma2[j]
  # }
  
  sigma2 <- rep(NaN, p)
  for(j in 1:p){
    sigma2[j] = object$fit[[j]]$sigma2 
  }
  var.error.hat <- object$var.error.hat
  var.error.new <- object$var.error.new
  
  res <- list(fit = object$fit, Theta = Theta_all, 
              Theta0 = Theta0, Theta1 = Theta1,
              Theta_U = Theta_U, Theta_W = Theta_W,
              Chat.new = Chat, Cor.new = Cor, 
              npc = npc, eigenfunctions = eigenfunctions, 
              eigenvalues = eigenvalues, U = Eig$vectors[,1:npc],
              npc0 = npc0, eigenfunctions0 = eigenfunctions0, 
              eigenvalues0 = eigenvalues0, U0 = Eig0$vectors[,1:npc0],
              npc1 = npc1, eigenfunctions1 = eigenfunctions1, 
              eigenvalues1 = eigenvalues1, U1 = Eig1$vectors[,1:npc1],
              var.error.hat = var.error.hat, 
              var.error.new = var.error.new, pve=pve, 
              G_invhalf = object$G_invhalf, sigma2 = sigma2, Beta = Beta)
  class(res) <- "mface.sparse"
  return(res)
  
}


sFPCA_post_reduced <- function(S0_fit, S1_fit, 
                               mfaces_fit, Beta, pve = 0.99){
  
  p <- 5
  
  S0 = S0_fit$S
  S1 = S1_fit$S
  
  beta1 = Beta[1]
  beta2 = Beta[2]
  beta3 = Beta[3]
  beta4 = Beta[4]
  beta5 = Beta[5]
  
  fit1 = mfaces_fit$fit$y1
  fit2 = mfaces_fit$fit$y2
  fit3 = mfaces_fit$fit$y3
  fit4 = mfaces_fit$fit$y4
  fit5 = mfaces_fit$fit$y5
  
  object = mfaces_fit
  tnew = object$argvals.new
  Bnew <- vector("list", p)
  B <- vector("list", p)
  for(k in 1:p){
    Bnew[[k]] <- object$fit[[k]]$Bnew %*% object$G_invhalf
    B[[k]] <- object$fit[[k]]$B %*% object$G_invhalf
  }
  Bbasis = object$fit[[k]]$Bnew %*% object$G_invhalf
  Bn = Bnew[[1]] / nrow(object$var.error.new)
  
  Theta0 = t(Bn) %*% S0 %*% Bn
  Theta1 = t(Bn) %*% S1 %*% Bn
  # Theta2 = t(Bn) %*% S2 %*% Bn
  # Theta3 = t(Bn) %*% S3 %*% Bn
  # Theta4 = t(Bn) %*% S4 %*% Bn
  # Theta5 = t(Bn) %*% S5 %*% Bn
  
  Theta11 = beta1^2*(Theta0 + Theta1)
  Theta12 = beta1*beta2*Theta0
  Theta13 = beta1*beta3*Theta0
  Theta14 = beta1*beta4*Theta0
  Theta15 = beta1*beta5*Theta0
  
  Theta21 = beta2*beta1*Theta0
  Theta22 = beta2^2*Theta0
  Theta23 = beta2*beta3*Theta0
  Theta24 = beta2*beta4*Theta0
  Theta25 = beta2*beta5*Theta0
  
  Theta31 = beta3*beta1*Theta0
  Theta32 = beta3*beta2*Theta0
  Theta33 = beta3^2*Theta0
  Theta34 = beta3*beta4*Theta0
  Theta35 = beta3*beta5*Theta0
  
  Theta41 = beta4*beta1*Theta0
  Theta42 = beta4*beta2*Theta0
  Theta43 = beta4*beta3*Theta0
  Theta44 = beta4^2*Theta0
  Theta45 = beta4*beta5*Theta0
  
  Theta51 = beta5*beta1*Theta0
  Theta52 = beta5*beta2*Theta0
  Theta53 = beta5*beta3*Theta0
  Theta54 = beta5*beta4*Theta0
  Theta55 = beta5^2*(Theta0 + Theta1)
  
  
  Theta_all = rbind(cbind(Theta11, Theta12, Theta13, Theta14, Theta15),
                    cbind(Theta21, Theta22, Theta23, Theta24, Theta25),
                    cbind(Theta31, Theta32, Theta33, Theta34, Theta35),
                    cbind(Theta41, Theta42, Theta43, Theta44, Theta45),
                    cbind(Theta51, Theta52, Theta53, Theta54, Theta55))
  Bnew <- do.call(bdiag, Bnew)#%*%G_invhalfmat
  Chat <- Bnew%*%Theta_all%*%t(Bnew)
  
  Theta_W = bdiag(beta1^2*Theta1, beta2^2*Theta1, beta3^2*Theta1, beta4^2*Theta1, beta5^2*Theta1)
  Theta_U = Theta_all - Theta_W
  
  
  Chat_diag = as.vector(diag(Chat))  
  Cor = diag(1/sqrt(Chat_diag))%*%Chat%*%diag(1/sqrt(Chat_diag))
  
  B <- do.call(bdiag, B)
  Chat.diag.pred = diag(as.matrix(tcrossprod(B%*%Matrix(Theta_all),B)))
  
  Eig <- eigen(Theta_all)
  npc <- which.max(cumsum(Eig$values)/sum(Eig$values)>pve)[1]
  
  eigenfunctions = matrix(Bnew%*%Eig$vectors[,1:min(npc, p*length(tnew))],
                          ncol=min(npc, p*length(tnew)))
  eigenvalues = Eig$values[1:min(npc, p*length(tnew))]
  
  Eig0 <- eigen(Theta0)
  npc0 <- which.max(cumsum(Eig0$values)/sum(Eig0$values)>pve)[1]
  
  eigenfunctions0 = matrix(Bbasis%*%Eig0$vectors[,1:min(npc0, length(tnew))],
                           ncol=min(npc0, length(tnew)))
  eigenvalues0 = Eig0$values[1:min(npc0, length(tnew))]
  
  Eig1 <- eigen(Theta1)
  npc1 <- which.max(cumsum(Eig1$values)/sum(Eig1$values)>pve)[1]
  
  eigenfunctions1 = matrix(Bbasis%*%Eig1$vectors[,1:min(npc1, length(tnew))],
                           ncol=min(npc1, length(tnew)))
  eigenvalues1 = Eig1$values[1:min(npc1, length(tnew))]
  
  # r_j = matrix(NaN, nrow(data),p)
  # for(j in 1:p){
  #   obs_j = data[,-c(1,2)][,j]
  #   r_j[,j] = (obs_j - object$fit[[j]]$mu.hat)^2
  # }
  # 
  # sigma2 = colMeans(matrix(c(r_j) - Chat.diag.pred, nrow(data), p))
  # sigma2[sigma2 <= 0] <- 1e-6
  # 
  # var.error.hat <- t(matrix(rep(sigma2, p*nrow(data)),p,nrow(data)))
  # var.error.new <- t(matrix(rep(sigma2, p*length(tnew)),p,length(tnew)))
  # 
  # for(j in 1:p){
  #   object$fit[[j]]$sigma2 = sigma2[j]
  # }
  
  sigma2 <- rep(NaN, p)
  for(j in 1:p){
    sigma2[j] = object$fit[[j]]$sigma2 
  }
  var.error.hat <- object$var.error.hat
  var.error.new <- object$var.error.new
  
  res <- list(fit = object$fit, Theta = Theta_all, 
              Theta0 = Theta0, Theta1 = Theta1,
              Theta_U = Theta_U, Theta_W = Theta_W,
              Chat.new = Chat, Cor.new = Cor, 
              npc = npc, eigenfunctions = eigenfunctions, 
              eigenvalues = eigenvalues, U = Eig$vectors[,1:npc],
              npc0 = npc0, eigenfunctions0 = eigenfunctions0, 
              eigenvalues0 = eigenvalues0, U0 = Eig0$vectors[,1:npc0],
              npc1 = npc1, eigenfunctions1 = eigenfunctions1, 
              eigenvalues1 = eigenvalues1, U1 = Eig1$vectors[,1:npc1],
              var.error.hat = var.error.hat, 
              var.error.new = var.error.new, pve=pve, 
              G_invhalf = object$G_invhalf, sigma2 = sigma2, Beta = Beta)
  class(res) <- "mface.sparse"
  return(res)
  
}



sFPCA_post_reduced_2 <- function(S0_fit, S1_fit, 
                               mfaces_fit, Beta, pve = 0.99){
  
  p <- 5
  
  S0 = S0_fit$S
  S1 = S1_fit$S
  
  beta1 = Beta[1]
  beta2 = Beta[2]
  beta3 = Beta[3]
  beta4 = Beta[4]
  beta5 = Beta[5]
  
  fit1 = mfaces_fit$fit$y1
  fit2 = mfaces_fit$fit$y2
  fit3 = mfaces_fit$fit$y3
  fit4 = mfaces_fit$fit$y4
  fit5 = mfaces_fit$fit$y5
  
  object = mfaces_fit
  tnew = object$argvals.new
  Bnew <- vector("list", p)
  B <- vector("list", p)
  for(k in 1:p){
    Bnew[[k]] <- object$fit[[k]]$Bnew %*% object$G_invhalf
    B[[k]] <- object$fit[[k]]$B %*% object$G_invhalf
  }
  Bbasis = object$fit[[k]]$Bnew %*% object$G_invhalf
  Bn = Bnew[[1]] / nrow(object$var.error.new)
  
  Theta0 = t(Bn) %*% S0 %*% Bn
  Theta1 = t(Bn) %*% S1 %*% Bn
  # Theta2 = t(Bn) %*% S2 %*% Bn
  # Theta3 = t(Bn) %*% S3 %*% Bn
  # Theta4 = t(Bn) %*% S4 %*% Bn
  # Theta5 = t(Bn) %*% S5 %*% Bn
  
  Theta11 = beta1^2*Theta0 
  Theta12 = beta1*beta2*Theta0
  Theta13 = beta1*beta3*Theta0
  Theta14 = beta1*beta4*Theta0
  Theta15 = beta1*beta5*Theta0
  
  Theta21 = beta2*beta1*Theta0
  Theta22 = beta2^2*Theta0
  Theta23 = beta2*beta3*Theta0
  Theta24 = beta2*beta4*Theta0
  Theta25 = beta2*beta5*Theta0
  
  Theta31 = beta3*beta1*Theta0
  Theta32 = beta3*beta2*Theta0
  Theta33 = beta3^2*Theta0
  Theta34 = beta3*beta4*Theta0
  Theta35 = beta3*beta5*Theta0
  
  Theta41 = beta4*beta1*Theta0
  Theta42 = beta4*beta2*Theta0
  Theta43 = beta4*beta3*Theta0
  Theta44 = beta4^2*Theta0
  Theta45 = beta4*beta5*Theta0
  
  Theta51 = beta5*beta1*Theta0
  Theta52 = beta5*beta2*Theta0
  Theta53 = beta5*beta3*Theta0
  Theta54 = beta5*beta4*Theta0
  Theta55 = beta5^2*(Theta0 + Theta1)
  
  
  Theta_all = rbind(cbind(Theta11, Theta12, Theta13, Theta14, Theta15),
                    cbind(Theta21, Theta22, Theta23, Theta24, Theta25),
                    cbind(Theta31, Theta32, Theta33, Theta34, Theta35),
                    cbind(Theta41, Theta42, Theta43, Theta44, Theta45),
                    cbind(Theta51, Theta52, Theta53, Theta54, Theta55))
  Bnew <- do.call(bdiag, Bnew)#%*%G_invhalfmat
  Chat <- Bnew%*%Theta_all%*%t(Bnew)
  
  Theta_W = bdiag(beta1^2*Theta1, beta2^2*Theta1, beta3^2*Theta1, beta4^2*Theta1, beta5^2*Theta1)
  Theta_U = Theta_all - Theta_W
  
  
  Chat_diag = as.vector(diag(Chat))  
  Cor = diag(1/sqrt(Chat_diag))%*%Chat%*%diag(1/sqrt(Chat_diag))
  
  B <- do.call(bdiag, B)
  Chat.diag.pred = diag(as.matrix(tcrossprod(B%*%Matrix(Theta_all),B)))
  
  Eig <- eigen(Theta_all)
  npc <- which.max(cumsum(Eig$values)/sum(Eig$values)>pve)[1]
  
  eigenfunctions = matrix(Bnew%*%Eig$vectors[,1:min(npc, p*length(tnew))],
                          ncol=min(npc, p*length(tnew)))
  eigenvalues = Eig$values[1:min(npc, p*length(tnew))]
  
  Eig0 <- eigen(Theta0)
  npc0 <- which.max(cumsum(Eig0$values)/sum(Eig0$values)>pve)[1]
  
  eigenfunctions0 = matrix(Bbasis%*%Eig0$vectors[,1:min(npc0, length(tnew))],
                           ncol=min(npc0, length(tnew)))
  eigenvalues0 = Eig0$values[1:min(npc0, length(tnew))]
  
  Eig1 <- eigen(Theta1)
  npc1 <- which.max(cumsum(Eig1$values)/sum(Eig1$values)>pve)[1]
  
  eigenfunctions1 = matrix(Bbasis%*%Eig1$vectors[,1:min(npc1, length(tnew))],
                           ncol=min(npc1, length(tnew)))
  eigenvalues1 = Eig1$values[1:min(npc1, length(tnew))]
  
  # r_j = matrix(NaN, nrow(data),p)
  # for(j in 1:p){
  #   obs_j = data[,-c(1,2)][,j]
  #   r_j[,j] = (obs_j - object$fit[[j]]$mu.hat)^2
  # }
  # 
  # sigma2 = colMeans(matrix(c(r_j) - Chat.diag.pred, nrow(data), p))
  # sigma2[sigma2 <= 0] <- 1e-6
  # 
  # var.error.hat <- t(matrix(rep(sigma2, p*nrow(data)),p,nrow(data)))
  # var.error.new <- t(matrix(rep(sigma2, p*length(tnew)),p,length(tnew)))
  # 
  # for(j in 1:p){
  #   object$fit[[j]]$sigma2 = sigma2[j]
  # }
  
  sigma2 <- rep(NaN, p)
  for(j in 1:p){
    sigma2[j] = object$fit[[j]]$sigma2 
  }
  var.error.hat <- object$var.error.hat
  var.error.new <- object$var.error.new
  
  res <- list(fit = object$fit, Theta = Theta_all, 
              Theta0 = Theta0, Theta1 = Theta1,
              Theta_U = Theta_U, Theta_W = Theta_W,
              Chat.new = Chat, Cor.new = Cor, 
              npc = npc, eigenfunctions = eigenfunctions, 
              eigenvalues = eigenvalues, U = Eig$vectors[,1:npc],
              npc0 = npc0, eigenfunctions0 = eigenfunctions0, 
              eigenvalues0 = eigenvalues0, U0 = Eig0$vectors[,1:npc0],
              npc1 = npc1, eigenfunctions1 = eigenfunctions1, 
              eigenvalues1 = eigenvalues1, U1 = Eig1$vectors[,1:npc1],
              var.error.hat = var.error.hat, 
              var.error.new = var.error.new, pve=pve, 
              G_invhalf = object$G_invhalf, sigma2 = sigma2, Beta = Beta)
  class(res) <- "mface.sparse"
  return(res)
  
}


sFPCA_post_all_S0 <- function(S0_fit,
                              mfaces_fit, Beta, pve = 0.99){
  
  p <- 5
  
  S0 = S0_fit$S
  
  beta1 = Beta[1]
  beta2 = Beta[2]
  beta3 = Beta[3]
  beta4 = Beta[4]
  beta5 = Beta[5]
  
  fit1 = mfaces_fit$fit$y1
  fit2 = mfaces_fit$fit$y2
  fit3 = mfaces_fit$fit$y3
  fit4 = mfaces_fit$fit$y4
  fit5 = mfaces_fit$fit$y5
  
  object = mfaces_fit
  tnew = object$argvals.new
  Bnew <- vector("list", p)
  B <- vector("list", p)
  for(k in 1:p){
    Bnew[[k]] <- object$fit[[k]]$Bnew %*% object$G_invhalf
    B[[k]] <- object$fit[[k]]$B %*% object$G_invhalf
  }
  Bbasis = object$fit[[k]]$Bnew %*% object$G_invhalf
  Bn = Bnew[[1]] / nrow(object$var.error.new)
  
  Theta0 = t(Bn) %*% S0 %*% Bn
  # Theta2 = t(Bn) %*% S2 %*% Bn
  # Theta3 = t(Bn) %*% S3 %*% Bn
  # Theta4 = t(Bn) %*% S4 %*% Bn
  # Theta5 = t(Bn) %*% S5 %*% Bn
  
  Theta11 = beta1^2*Theta0
  Theta12 = beta1*beta2*Theta0
  Theta13 = beta1*beta3*Theta0
  Theta14 = beta1*beta4*Theta0
  Theta15 = beta1*beta5*Theta0
  
  Theta21 = beta2*beta1*Theta0
  Theta22 = beta2^2*Theta0
  Theta23 = beta2*beta3*Theta0
  Theta24 = beta2*beta4*Theta0
  Theta25 = beta2*beta5*Theta0
  
  Theta31 = beta3*beta1*Theta0
  Theta32 = beta3*beta2*Theta0
  Theta33 = beta3^2*Theta0
  Theta34 = beta3*beta4*Theta0
  Theta35 = beta3*beta5*Theta0
  
  Theta41 = beta4*beta1*Theta0
  Theta42 = beta4*beta2*Theta0
  Theta43 = beta4*beta3*Theta0
  Theta44 = beta4^2*Theta0
  Theta45 = beta4*beta5*Theta0
  
  Theta51 = beta5*beta1*Theta0
  Theta52 = beta5*beta2*Theta0
  Theta53 = beta5*beta3*Theta0
  Theta54 = beta5*beta4*Theta0
  Theta55 = beta5^2*Theta0
  
  
  Theta_all = rbind(cbind(Theta11, Theta12, Theta13, Theta14, Theta15),
                    cbind(Theta21, Theta22, Theta23, Theta24, Theta25),
                    cbind(Theta31, Theta32, Theta33, Theta34, Theta35),
                    cbind(Theta41, Theta42, Theta43, Theta44, Theta45),
                    cbind(Theta51, Theta52, Theta53, Theta54, Theta55))
  Bnew <- do.call(bdiag, Bnew)#%*%G_invhalfmat
  Chat <- Bnew%*%Theta_all%*%t(Bnew)
  
  Chat_diag = as.vector(diag(Chat))  
  Cor = diag(1/sqrt(Chat_diag))%*%Chat%*%diag(1/sqrt(Chat_diag))
  
  B <- do.call(bdiag, B)
  Chat.diag.pred = diag(as.matrix(tcrossprod(B%*%Matrix(Theta_all),B)))
  
  Eig <- eigen(Theta_all)
  npc <- which.max(cumsum(Eig$values)/sum(Eig$values)>pve)[1]
  
  eigenfunctions = matrix(Bnew%*%Eig$vectors[,1:min(npc, p*length(tnew))],
                          ncol=min(npc, p*length(tnew)))
  eigenvalues = Eig$values[1:min(npc, p*length(tnew))]
  
  Eig0 <- eigen(Theta0)
  npc0 <- which.max(cumsum(Eig0$values)/sum(Eig0$values)>pve)[1]
  
  eigenfunctions0 = matrix(Bbasis%*%Eig0$vectors[,1:min(npc0, length(tnew))],
                           ncol=min(npc0, length(tnew)))
  eigenvalues0 = Eig0$values[1:min(npc0, length(tnew))]
  # r_j = matrix(NaN, nrow(data),p)
  # for(j in 1:p){
  #   obs_j = data[,-c(1,2)][,j]
  #   r_j[,j] = (obs_j - object$fit[[j]]$mu.hat)^2
  # }
  # 
  # sigma2 = colMeans(matrix(c(r_j) - Chat.diag.pred, nrow(data), p))
  # sigma2[sigma2 <= 0] <- 1e-6
  # 
  # var.error.hat <- t(matrix(rep(sigma2, p*nrow(data)),p,nrow(data)))
  # var.error.new <- t(matrix(rep(sigma2, p*length(tnew)),p,length(tnew)))
  # 
  # for(j in 1:p){
  #   object$fit[[j]]$sigma2 = sigma2[j]
  # }
  
  sigma2 <- rep(NaN, p)
  for(j in 1:p){
    sigma2[j] = object$fit[[j]]$sigma2 
  }
  var.error.hat <- object$var.error.hat
  var.error.new <- object$var.error.new
  
  res <- list(fit = object$fit, Theta = Theta_all, 
              Theta0 = Theta0,
              Chat.new = Chat, Cor.new = Cor, 
              npc = npc, eigenfunctions = eigenfunctions, 
              eigenvalues = eigenvalues, U = Eig$vectors[,1:npc],
              npc0 = npc0, eigenfunctions0 = eigenfunctions0, 
              eigenvalues0 = eigenvalues0, U0 = Eig0$vectors[,1:npc0],
              var.error.hat = var.error.hat, 
              var.error.new = var.error.new, pve=pve, 
              G_invhalf = object$G_invhalf, sigma2 = sigma2, Beta = Beta)
  class(res) <- "mface.sparse"
  return(res)
}
