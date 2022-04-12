## Plot the estimated mean and eigenfunctions
library(rstan)
library(tidyverse)
library(gridExtra)
library(latex2exp)
library(xtable)
library(orthogonalsplinebasis)

source('source_code/functions.R')
load('RData/summary/summary_posterior_M2_whole.RData')

f = function(dat, lab.x, lab.y, range.y){
  p <- ggplot(data = dat, aes(x = x, y = y)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'grey70') + 
    geom_line(aes(y = y), linetype = 'solid', size = 1) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill=NA, size=1),
          plot.margin = margin(0.5, 0.25, 0, 0.25, "cm")) +
    xlab(TeX(lab.x)) +
    ylab(TeX(lab.y)) + 
    ylim(range.y[1], range.y[2])
  return(p)
}

f2 = function(dat, lab.x, lab.y, range.y){
  p <- ggplot(data = dat, aes(x = x, y = y)) + 
    geom_line(linetype = 'solid', size = 1) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill=NA, size=1),
          plot.margin = margin(0.5, 0.25, 0, 0.25, "cm")) +
    xlab(TeX(lab.x)) +
    ylab(TeX(lab.y)) + 
    ylim(range.y[1], range.y[2])
  return(p)
}

J = 5
P = 9
tnew = (0:101)/10

summ = read.csv("RData/summary/summary_posterior_M2_whole.csv")

bs.mat = create_bs(time.grid = tnew, pred.time = tnew, nbasis = P)

range.y <- list(y1 = c(0, 15), y2 = c(0, 10), y3 = c(0, 10), 
                y4 = c(0, 5), y5 = c(0, 5))
mu <- vector('list', length = J)
p = list()
for (j in 1:J){
  B_mean = summ[((j-1)*P+1):(j*P), 2]
  B_lower = summ[((j-1)*P+1):(j*P), 4]
  B_upper = summ[((j-1)*P+1):(j*P), 5]
  B_est = (bs.mat %*% B_mean)[, 1]
  mu[[j]] = B_est
  dat = data.frame(x = tnew, y = B_est, lower = bs.mat %*% B_lower, upper = bs.mat %*% B_upper)
  lab.y = paste0('$\\mu_', j, '(t)$')
  p[[j]] <- f(dat, 'Time', lab.y, range.y[[j]])
}

# Plot eigenfunctions ##
dat = data.frame(x = tnew, y = phi_est[, 1])
p[[6]] = f2(dat, 'Time', '$\\phi_1(t)$', c(0.5, 1.5))
dat = data.frame(x = tnew, y = phi_est[, 2])
p[[7]] = f2(dat, 'Time', '$\\phi_2(t)$', c(-2.5, 2))
dat = data.frame(x = tnew, y = psi_est[, 1])
p[[8]] = f2(dat, 'Time', '$\\psi_1(t)$', c(-0.5, 3))

setEPS()
postscript("plot/est_func.eps", height = 8, width = 12)
grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], p[[7]], p[[8]], 
             nrow = 2)
dev.off()

save(list = c('mu', 'phi_est', 'psi_est'), file = 'RData/true_function.RData')
