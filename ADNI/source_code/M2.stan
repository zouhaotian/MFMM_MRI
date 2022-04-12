data{
  int n; // number of subjects
  int J; // number of longitudinal outcomes
  int nobs; // number of observations
  int nmiss[J]; // number of missingness
  int id_long[nobs]; // ID
  int L0; // number of eigenvalues for U_i(t)
  int L1; // number of eigenvalues for W_ij(t)
  int Lm; // number of eigenfunctions for f_mi(v)
  int P; // number of cubic B-spline functions
  int P_surv; 
  
  real Y1[nobs]; // Y1
  real Y2[nobs]; // Y2
  real Y3[nobs];
  real Y4[nobs];
  real Y5[nobs];
  int miss_index1[nmiss[1]]; 
  int miss_index2[nmiss[2]]; 
  int miss_index3[nmiss[3]]; 
  int miss_index4[nmiss[4]]; 
  int miss_index5[nmiss[5]]; 
  real time[nobs]; // observed time
  
  matrix[n, P_surv] x; 
  
  real surv_time[n]; // survival time
  real status[n]; // censoring status
  matrix[nobs, P] b; // cubic B-spline matrix
  matrix[nobs, L0] phi; // estimated eigenfunctions for U_i(t)
  matrix[nobs, L1] psi; // estimated eigenfunctions for W_ij(t)
  matrix[n, Lm] m_mat; // m_{ij} matrix
  matrix[L0, Lm] f_l; // \int_V phi_{ml}(v) psi_{mj}(v) dv
  real beta_m; 
  matrix[n, L0] phi_surv; // phi_l(T_i) 
  matrix[n, L1] psi_surv; // psi_l(T_i)
  
  int<lower=0> Ltau; // number of knots for baseline hazard function
  real tau[Ltau]; // knots for time 
  matrix[n, Ltau+1] h_grid; // observed survival time covariate for survival outcome
  matrix[n, Ltau] h_index; // observed time spline coefficient index for survival outcome
  
  int<lower=0> G; // number of quadrature points
  vector[G] w; // weights for quadrature
  matrix[n, Ltau] constant; // (b[i, e] - a[i, e])/2
  row_vector[G] phi1_interval[n, Ltau]; // phi_1 evaluated at G quadrature points for E intervals
  row_vector[G] phi2_interval[n, Ltau]; // phi_2 evaluated at G quadrature points for E intervals
  row_vector[G] psi1_interval[n, Ltau]; // psi_1 evaluated at G quadrature points for E intervals
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
  vector<lower=0>[J] sigma; 
  
  vector[Ltau] logh0; 
  vector[P_surv] gamma_x; // coefficients for x
  real gamma0; // coefficients for U_i(t)
  real gamma1[J];  // coefficients for W_{ij}(t)
  vector[Lm] gamma_m; // coefficients for zeta_{mil}
  real Y1_imp[nmiss[1]]; 
  real Y2_imp[nmiss[2]];
  real Y3_imp[nmiss[3]];
  real Y4_imp[nmiss[4]]; 
  real Y5_imp[nmiss[5]];
}
transformed parameters{
  real Y1_full[nobs] = Y1; 
  real Y2_full[nobs] = Y2;
  real Y3_full[nobs] = Y3;
  real Y4_full[nobs] = Y4;
  real Y5_full[nobs] = Y5;
  
  real mu1[nobs];
  real mu2[nobs];
  real mu3[nobs];
  real mu4[nobs];
  real mu5[nobs];
  matrix[n, Lm] mu_matrix = beta_m*(xi'*f_l); // beta_m * (i,j-th element is sum_l xi_{il}*f_{lj})
  matrix[n, Lm] xi_m = m_mat - mu_matrix; // FPC scores for f_mi(v)
  
  vector[L0] d0;
  vector[L1] d1;
  vector[L0] tmp_xi;
  
  vector[L1] tmp_zeta1;
  vector[L1] tmp_zeta2;
  vector[L1] tmp_zeta3;
  vector[L1] tmp_zeta4;
  vector[L1] tmp_zeta5;
  
  real rs_w[n]; // risk score for subject i
  real h[n]; // hazard function
  matrix[n, Ltau] H; // cumulative hazard in interval (e, e+1)
  real LL[n]; // log survival likelihood
  row_vector[G] U; // for subject i, interval e, U = tmp_xi[1]*phi_surv[i, e] 
  row_vector[G] W1; // for subject i, interval e, W1 = tmp_zeta1*psi_surv[i, e]
  row_vector[G] W2;
  row_vector[G] W3;
  row_vector[G] W4;
  row_vector[G] W5;
  row_vector[G] f; // f = exp(\gamma_0*U_i(t) + \sum_{j=1}^J \gamma_{1j}*W_{ij}(t))
  
  Y1_full[miss_index1] = Y1_imp;
  Y2_full[miss_index2] = Y2_imp;
  Y3_full[miss_index3] = Y3_imp;
  Y4_full[miss_index4] = Y4_imp;
  Y5_full[miss_index5] = Y5_imp;
  
  for (i in 1:L0){
    d0[i] = pow(sqrt_d0[i], 2);
  }
  
  for (i in 1:L1){
    d1[i] = pow(sqrt_d1[i], 2);
  }
  
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
  
  for (i in 1:n){
    tmp_xi = col(xi, i);
    
    tmp_zeta1 = col(zeta_1, i);
    tmp_zeta2 = col(zeta_2, i);
    tmp_zeta3 = col(zeta_3, i);
    tmp_zeta4 = col(zeta_4, i);
    tmp_zeta5 = col(zeta_5, i);
    
    rs_w[i] = gamma0*(phi_surv[i]*tmp_xi) + gamma1[1]*(psi_surv[i]*tmp_zeta1) + 
      gamma1[2]*(psi_surv[i]*tmp_zeta2) + gamma1[3]*(psi_surv[i]*tmp_zeta3) + 
      gamma1[4]*(psi_surv[i]*tmp_zeta4) + gamma1[5]*(psi_surv[i]*tmp_zeta5);
    h[i] = exp(h_index[i]*logh0 + x[i]*gamma_x + xi_m[i]*gamma_m + rs_w[i]);
    
    for (k in 1:Ltau){
      U = tmp_xi[1]*phi1_interval[i, k] + tmp_xi[2]*phi2_interval[i, k];
      W1 = tmp_zeta1[1]*psi1_interval[i, k];
      W2 = tmp_zeta2[1]*psi1_interval[i, k];
      W3 = tmp_zeta3[1]*psi1_interval[i, k];
      W4 = tmp_zeta4[1]*psi1_interval[i, k];
      W5 = tmp_zeta5[1]*psi1_interval[i, k];
      f = exp(gamma0*U + gamma1[1]*W1 + gamma1[2]*W2 + gamma1[3]*W3 + gamma1[4]*W4 + gamma1[5]*W5);
      H[i, k] = exp(x[i]*gamma_x + logh0[k] + xi_m[i]*gamma_m)*constant[i, k]*(f*w);
    }
    LL[i] = status[i]*log(h[i]) + sum(-H[i]);
  }
}
model{
  A1 ~ normal(0, 10);
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
  
  logh0 ~ normal(0, 10);
  gamma_x ~ normal(0, 10);
  gamma0 ~ normal(0, 10);
  gamma1 ~ normal(0, 10);
  gamma_m ~ normal(0, 10); 
  
  target+=skew_normal_lpdf(Y1_full | mu1, sigma[1], 0);
  target+=skew_normal_lpdf(Y2_full | mu2, sigma[2], 0);
  target+=skew_normal_lpdf(Y3_full | mu3, sigma[3], 0);
  target+=skew_normal_lpdf(Y4_full | mu4, sigma[4], 0);
  target+=skew_normal_lpdf(Y5_full | mu5, sigma[5], 0);
  target+=LL;
}
generated quantities{
  vector[Lm] dm; 
  for (i in 1:Lm){
    dm[i] = variance(col(xi_m, i)); 
  }
}
