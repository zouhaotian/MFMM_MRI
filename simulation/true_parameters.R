sim.model = c("s1_2S", "s1_bCox", "s1_MFMM_MRI_x2", "s1_MFMM_MRI_x3", "s1_MFMM_MRI_x4", 
              "s1_MJM_MRI", "s1_NM", "s2_MFMM_MRI", "s3_MFMM_MRI")

load('RData/data_gen1.RData')

tp$d0 = c(2, 1)
tp$d1 = 1
tp$beta_m = -0.05
tp$omega_m = 0.4
tp$d_m = c(0.10, 0.07, 0.05, 0.03, 0.01)

save(tp, file = 'RData/data_genx.RData')

load('RData/data_genx.RData')
tp$N.train = 1200
tp$N.test = 400
save(tp, file = 'RData/data_genx2.RData')


load('RData/data_genx.RData')
tp$logh0 = -1
save(tp, file = 'RData/data_genx3.RData')


load('RData/data_genx.RData')
tp$N.train = 1200
tp$N.test = 400
tp$logh0 = -1
save(tp, file = 'RData/data_genx4.RData')

load('RData/data_genx.RData')
tp$beta_m = 0
save(tp, file = 'RData/data_gen2.RData')


load('RData/data_genx.RData')
save(tp, file = 'RData/data_gen3.RData')
