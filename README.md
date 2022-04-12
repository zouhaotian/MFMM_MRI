# MFMM_MRI
This is the code repository for multivariate functional mixed model with MRI data: an application to Alzheimer’s Disease.

There are two folders: ADNI for real data application; simulation for the simulation study.

In the ADNI folder, please follow the steps to reproduce results:

(1) Run data_query.R to generate dataset/ADNI_long.csv and dataset/ADNI_surv.csv. The raw dataset are not placed in the dataset folder due to privacy.

(2) Run voxel_selection.R and get_voxel.R to generate dataset/Surv.mri.dat.csv. Run get_hippocampus_voxel.R to generate dataset/hippo.mri.dat.csv (not placed due to excessive file size).

(3) Run FPCA_hippo.R to generate RData/FPCA_Hippo_result.RData. Run FPCA_whole.R to generate RData/FPCA_whole_result.RData.

(4) Run ADNI_data_transform.R to generate dataset/ADNI_long_2.csv, for Box-Cox transformation.

(5) Run ADNI_spaghetti.R to generate Figure 1 (plot/spaghetti.eps).

(6) Run ADNI_M1_whole.R, ADNI_M1.R, ADNI_M1_MFMM.R to fit models M1-whole, M1-hippo, and M1. We run similar models for M2.

(7) Run summary_MP.R to generate result/summary_MP.csv as Table 1.

(8) Run brain_heat_map.R and use BrainNet Viewer to generate Figure 3 (plot/brain_heat_map_c3.eps). Run LRT_gamma_m.R to generate LRT_gamma_m.RData and LRT_gamma_m_Hippo.RData. Run summary_beta_m.R to generate RData/summary.beta_m.RData. Run summary_function.R to generate Figure 1 in the Supplementary Material (plot/est_func.eps).

(9) Run ADNI_M1_whole_dp.R, ADNI_M1_dp.R, ADNI_M1_MFMM_dp.R to perform internal cross-validation of dynamic prediction for Model 1. Similarly, we perform dynamic prediction for Model 2.

(10) Run summary_dp.R to generate the iAUC and iBS values (result/iAUC_iBS.csv).

(11) Run dp_plot.R to generate Figure 4 (plot/727_1.eps and plot/727_2.eps).


For the simulation study, please follow the steps to reproduce results:

(1) Run simulation1.R to perform 100 simulation replications.

(2) Run summary1_posterior.R to generate RData/summary1.csv and RData/integrated_AUC_BS1.csv.
