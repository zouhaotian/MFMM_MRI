library(tidyverse)
library(rstan)
library(data.table)
library(EveTemplate)
library(orthogonalsplinebasis)

source('source_code/f1.R')

MRI.dat <- fread("dataset/Surv.mri.dat.csv") %>% as.matrix()
N = nrow(MRI.dat)
V = ncol(MRI.dat)
L0 = 2
Lm = 5
obsgrid = (1:V)/V

M.train = MRI.dat

Pm = min(round(V/30), 35)
mu_m = colMeans(M.train)
mu_m_smooth = bs.smooth(mu_m, obsgrid, obsgrid, nbasis = Pm)
mu_m_expand = matrix(rep(mu_m_smooth$est.value, N), nrow = N, byrow = T)
M.demean = M.train - mu_m_expand

fname <- paste0('RData/summary/summary_posterior_M2_whole.RData')
load(fname)

## Create m_{ij} matrix ##
m_mat = matrix(NA, nrow = N, ncol = Lm)
f_l = matrix(NA, nrow = L0, ncol = Lm)

for (j in 1:Lm){
  for (i in 1:N){
    m_mat[i, j] <- c_int(M.demean[i, ], Psi_m_est[, j], width = 1/(V-1))
  }
  for (l in 1:L0){
    f_l[l, j] <- c_int(Phi_m_est[, l], Psi_m_est[, j], width = 1/(V-1))
  }
}

fname <- paste0('stan_fit/M2_whole_posterior.RData')
load(fname)

xi_sample = extract(fitStan, pars = 'xi')$xi
xi <- colMeans(xi_sample)

## Express m_i(v) = mu_m(v) + beta_m*u_{mi}(v) + f_{mi}(v) + e_{mi}(v)
## Let i = 1
i <- 1

mi <- MRI.dat[i, ]
mu_m <- mu_m_smooth$est.value
u_mi <- (Phi_m_est %*% xi[, i])[, 1]

xi_mi <- m_mat[1, ] - beta_m* (t(f_l) %*% xi[, i])
f_mi <- (Psi_m_est %*% xi_mi)[, 1]

fname = 'RData/brain_heat_map.RData'
save(list = c('mi', 'mu_m', 'u_mi', 'f_mi'), file = fname)
# load(fname)

## rescale the eve template coordinates ##
fname <- 'dataset/Eve_coord.txt'
Eve.coord <- read.table(fname, sep = ' ', header = T)
Eve.coord <- Eve.coord %>% arrange(x, y, z)

eve_brain <- readEve('Brain')
eve_brain_mask <- readEve('Brain_Mask')
dat.mask <- img_data(eve_brain_mask)
dims <- prod(dim(dat.mask))
x.coord <- rep(1:181, each = 217*181)
y.coord <- rep(rep(1:217, each = 181), 181)
z.coord <- rep(1:181, 181*217)
dat2 <- data.table(x = x.coord, y = y.coord, z = z.coord,
                   index = c(1:dims), value = as.vector(dat.mask))
dat3 <- dat2[value==1, ] %>% select(-value)
## records mask index that corresponds to the whole brain index
mask.dat <- data.frame(index = dat2[eve_brain_mask==1, ]$index, 
                       mask.index = c(1:1843303)) 

## rescale the mask to range of ICBM152 surface ##
x.new.range <- range(Eve.coord$x)
y.new.range <- range(Eve.coord$y)
z.new.range <- range(Eve.coord$z)

rescale <- function(x, new.range){
  x.range <- range(x)
  x.diff <- x.range[2] - x.range[1]
  new.range.diff <- new.range[2] - new.range[1]
  x.new <- (x - x.range[1])*(new.range.diff/x.diff) + new.range[1]
  return(x.new)
}

x.new <- rescale(dat3$x, x.new.range)
y.new <- rescale(dat3$y, y.new.range)
z.new <- rescale(dat3$z, z.new.range)

label_df <- readEveMap('II')
eve_labels <- getEveMapLabels('II')
label_dat <- data.frame(index = c(1:dims), integer_label = as.vector(img_data(label_df))) %>% 
  left_join(eve_labels, by = 'integer_label') %>% select(index, integer_label, text_label)

eve.adjusted <- data.frame(x = x.new, y = y.new, z = z.new, 
                           index = dat3$index) %>% 
  left_join(label_dat, by = 'index')

load('RData/p_value_index.RData')
Surv.index <- index
surv.dat <- data.frame(index = mask.dat[mask.dat$mask.index %in% Surv.index, ]$index) %>% 
  left_join(eve.adjusted, by = 'index')

surv.dat$mi = mi
surv.dat$mu_m = mu_m
surv.dat$u_mi = u_mi
surv.dat$f_mi = f_mi

fname = 'RData/surv_brain_heat_map.RData'
save(list = 'surv.dat', file = fname)
# load(fname)

selected.label <- sort(unique(surv.dat$integer_label))[-1]
exclude.label <- c(3, 4, 5, 6, 7, 23)
selected.label <- selected.label[which(!selected.label %in% exclude.label)]
surv.dat3 <- surv.dat[surv.dat$integer_label %in% selected.label, ]

get_node_size <- function(tmp.dat, node.size = 0.1){
  tmp.dat2 <- tmp.dat
  for (i in 1:nrow(tmp.dat2)){
    tmp.dat2$node.size[i] <- nrow(tmp.dat[tmp.dat$integer_label==tmp.dat2$integer_label[i], ])*node.size
  }
  return(tmp.dat2)
}

surv.dat4 <- get_node_size(surv.dat3, node.size = 0.2)
surv.dat4$nodup.label <- ifelse(duplicated(surv.dat4$integer_label), "-", surv.dat4$integer_label)

write.table(format(surv.dat4[c(2, 3, 4, 7, 11, 12)], digits = 2), 
            file = 'dataset/B_mi.node',
            row.names = F,
            col.names = F,
            quote = F)
write.table(format(surv.dat4[c(2, 3, 4, 8, 11, 12)], digits = 2), 
            file = 'dataset/B_mum.node',
            row.names = F,
            col.names = F,
            quote = F)
write.table(format(surv.dat4[c(2, 3, 4, 9, 11, 12)], digits = 2), 
            file = 'dataset/B_umi.node',
            row.names = F,
            col.names = F,
            quote = F)
write.table(format(surv.dat4[c(2, 3, 4, 10, 11, 12)], digits = 2), 
            file = 'dataset/B_fmi.node',
            row.names = F,
            col.names = F,
            quote = F)
