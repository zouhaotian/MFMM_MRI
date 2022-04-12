data{
  int n; // number of subjects
  int J; // number of longitudinal outcomes
  int nobs; // number of observations
  int nmiss[J]; // number of missingness
  int id_long[nobs]; // ID
  int L0; // number of eigenvalues for U_i(t)
  int L1; // number of eigenvalues for W_ij(t)
  int P; // number of cubic B-spline functions
  
  real Y1[nobs]; // Y1
  real Y2[nobs]; // Y2
  real Y3[nobs];
  real Y4[nobs];
  real Y5[nobs];
  int miss_index1[nmiss[1]];  // missing index for Y1
  int miss_index2[nmiss[2]];
  int miss_index3[nmiss[3]]; 
  int miss_index4[nmiss[4]];
  int miss_index5[nmiss[5]]; 
  real time[nobs]; // observed time

  matrix[nobs, P] b; // cubic B-spline matrix
  matrix[nobs, L0] phi; // estimated eigenfunctions for U_i(t)
  matrix[nobs, L1] psi; // estimated eigenfunctions for W_ij(t)
}
parameters{
  vector[P] A1; // coefficients for mu1
  vector[P] A2; // coefficients for mu2
  vector[P] A3;
  vector[P] A4;
  vector[P] A5;
  
  vector<lower=0>[L0] sqrt_d0; // standard deviation for FPC scores xi_il
  vector<lower=0>[L1] sqrt_d1; // standard deviation for zeta_ijl
  
  matrix[L0, n] xi; // FPC scores for U_i(t)
  matrix[L1, n] zeta_1; // FPC scores for W_i1(t)
  matrix[L1, n] zeta_2; // FPC scores for W_i2(t)
  matrix[L1, n] zeta_3;
  matrix[L1, n] zeta_4;
  matrix[L1, n] zeta_5;
  
  vector[J-1] beta; // beta
  vector<lower=0>[J] sigma; // sd of Y_{ij}(t)
  
  real Y1_imp[nmiss[1]]; 
  real Y2_imp[nmiss[2]];
  real Y3_imp[nmiss[3]];
  real Y4_imp[nmiss[4]];
  real Y5_imp[nmiss[5]];
}
transformed parameters{
  real Y1_full[nobs] = Y1; // Full Y1 (observed + imputed)
  real Y2_full[nobs] = Y2;
  real Y3_full[nobs] = Y3;
  real Y4_full[nobs] = Y4;
  real Y5_full[nobs] = Y5;
  
  real mu1[nobs]; // m_{i1}(t)
  real mu2[nobs];
  real mu3[nobs];
  real mu4[nobs];
  real mu5[nobs];
  
  vector[L0] d0 = sqrt_d0 .* sqrt_d0;
  vector[L1] d1 = sqrt_d1 .* sqrt_d1;
  vector[L0] tmp_xi;
  
  vector[L1] tmp_zeta1;
  vector[L1] tmp_zeta2;
  vector[L1] tmp_zeta3;
  vector[L1] tmp_zeta4;
  vector[L1] tmp_zeta5;

  Y1_full[miss_index1] = Y1_imp;
  Y2_full[miss_index2] = Y2_imp;
  Y3_full[miss_index3] = Y3_imp;
  Y4_full[miss_index4] = Y4_imp;
  Y5_full[miss_index5] = Y5_imp;
  
  for (i in 1:nobs){
    tmp_xi = col(xi, id_long[i]); // FPC scores xi_{i} for subject i
    
    tmp_zeta1 = col(zeta_1, id_long[i]); // FPC scores zeta_{i1} for subject i
    tmp_zeta2 = col(zeta_2, id_long[i]); // FPC scores zeta_{i2} for subject i
    tmp_zeta3 = col(zeta_3, id_long[i]);
    tmp_zeta4 = col(zeta_4, id_long[i]);
    tmp_zeta5 = col(zeta_5, id_long[i]);
    
    mu1[i] = b[i]*A1 + phi[i]*tmp_xi + psi[i]*tmp_zeta1;
    mu2[i] = b[i]*A2 + beta[1]*(phi[i]*tmp_xi + psi[i]*tmp_zeta2);
    mu3[i] = b[i]*A3 + beta[2]*(phi[i]*tmp_xi + psi[i]*tmp_zeta3);
    mu4[i] = b[i]*A4 + beta[3]*(phi[i]*tmp_xi + psi[i]*tmp_zeta4);
    mu5[i] = b[i]*A5 + beta[4]*(phi[i]*tmp_xi + psi[i]*tmp_zeta5);
  }
}
model{
  A1 ~ normal(0, 10); // prior 
  A2 ~ normal(0, 10);
  A3 ~ normal(0, 10);
  A4 ~ normal(0, 10);
  A5 ~ normal(0, 10);
  
  for (i in 1:L0){
    sqrt_d0[i] ~ inv_gamma(0.1, 0.1);
    xi[i] ~ normal(0, sqrt_d0[i]);
  }
  
  for (i in 1:L1){
    sqrt_d1[i] ~ inv_gamma(0.1, 0.1);
    zeta_1[i] ~ normal(0, sqrt_d1[i]);
    zeta_2[i] ~ normal(0, sqrt_d1[i]);
    zeta_3[i] ~ normal(0, sqrt_d1[i]);
    zeta_4[i] ~ normal(0, sqrt_d1[i]);
    zeta_5[i] ~ normal(0, sqrt_d1[i]);
  }
  
  beta ~ normal(0, 10);
  sigma ~ inv_gamma(0.1, 0.1);
  
  target+=skew_normal_lpdf(Y1_full | mu1, sigma[1], 0); // likelihood
  target+=skew_normal_lpdf(Y2_full | mu2, sigma[2], 0);
  target+=skew_normal_lpdf(Y3_full | mu3, sigma[3], 0);
  target+=skew_normal_lpdf(Y4_full | mu4, sigma[4], 0);
  target+=skew_normal_lpdf(Y5_full | mu5, sigma[5], 0);
}
