load('RData/summary/summary_posterior_M2_whole_MP.RData')
l_H1 = l$Dbar.surv/(-2)

load('RData/summary/summary_prior_M2_MP.RData')
l_H0 = l$Dbar.surv/(-2)

load('RData/FPCA_whole_result.RData')
Lm = length(FPCA.res$values) - 2 ## Lm = 5 for whole-brain voxels
l = -2*(l_H0 - l_H1) + Lm

## Under H0, l is asymptotically distributed with chi-sq with Lm df
p_value = 1 - pchisq(l, df = Lm)

save(list = c('l', 'l_H0', 'l_H1', 'Lm', 'p_value'), file = 'RData/LRT_gamma_m.RData')

## Test of hippo vs. no MRI
load('RData/summary/summary_posterior_M2_MP.RData')
l_H1 = l$Dbar.surv/(-2)

load('RData/summary/summary_prior_M2_MP.RData')
l_H0 = l$Dbar.surv/(-2)

load('RData/FPCA_Hippo_result.RData')
Lm = length(FPCA.res$values) - 2 
l = -2*(l_H0 - l_H1) + Lm

## Under H0, l is asymptotically distributed with chi-sq with Lm df
p_value = 1 - pchisq(l, df = Lm)

save(list = c('l', 'l_H0', 'l_H1', 'Lm', 'p_value'), file = 'RData/LRT_gamma_m_Hippo.RData')

